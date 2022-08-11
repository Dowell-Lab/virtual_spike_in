#!/usr/bin/env nextflow
/*
========================================================================================
                         VirtualSpikeInFlow - VSI Pipeline
========================================================================================
 Virtual Spike in Pipeline
 #### Homepage / Documentation
 https://github.com/Dowell-Lab/Analysis-Flow
 #### Authors
 Zachary Maas <zama8258@colorado.edu>
========================================================================================
========================================================================================
*/

def helpMessage() {
	log.info"""
	=========================================
	AnalysisFlow v${params.version}
  =========================================
		This pipeline requires that you manually create a configuration file to
		analyze your data, due to the large number of parameters involved.

	Usage:
  The typical command for running the pipeline is as follows:
		nextflow run main.nf -profile example
		Required arguments:
		-profile                      Configuration profile to use. <base, slurm, example>
		--workdir                     Nextflow working directory where all intermediate files are saved.
		""".stripIndent()
	}

params.help = false
if (params.help){
    helpMessage()
    exit 0
}

println "[Log]: Threshold is ${params.threshold}"

// Step 0 -- Convert cram to bam if needed
if ( params.cramDir ) {
		// We generate cram names from bedgraph names to ensure that they match.
		println "[Log]: Creating BAM files from CRAM files..."
		samples = Channel
			.fromPath(params.cramDir)
			.map { file -> tuple(file.baseName, file)}

		process cramToBam {
		cpus 16
		memory '4 GB'
		time '30m'
		tag "$prefix"

		input:
				set val(prefix), file(cram) from samples

		output:
				set val(prefix), file("${prefix}.sorted.bam") into bamSamples
				file("${prefix}.sorted.bam") into bamInit
				file("${prefix}.sorted.bam") into bamForCounts

				module 'samtools'
				script:
				"""
				samtools view -@ 16 -b -1 -T ${params.reffasta} ${cram} > ${prefix}.sorted.bam
				"""
		}
} else {
		// We generate bam names from bedgraph names to ensure that they match.
		// TODO Not necessary
		samples = Channel
		.fromPath(params.bedgraphs)
		.map { file -> tuple(file.baseName, file, "$params.bamDir" + "$file.baseName" + ".sorted.bam")}
		.set() { bamSamples }

		// Match what we do in the previous step in terms of variable names
		(bamInit, bamForCounts) = Channel
				.fromPath(params.bedgraphs)
				.map { it -> file("$params.bamDir" + "$it.baseName" + ".sorted.bam") }
				.into(2)
}

// Step 1 -- Filter for maximal isoforms
if (params.filter == "yes") {
		process filterIsoform {
				cpus 4
				memory '4GB'
				time '1h'
				errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
				maxRetries 3
				maxErrors -1
				tag "$prefix"
				publishDir "${params.outdir}/isoform/", mode: 'copy', pattern: "*.sorted.isoform_max.bed", overwrite: true
				input:
				set val(prefix), val(bam) from bamSamples

				output:
				set val(prefix), file("*.sorted.isoform_max.bed") into filteredIsoforms

				module 'python/3.6.3'
				module 'bedtools'
				module 'subread'
				script:
				"""
				export NUM_CORES=4
				export InterestFile=${bam}
				export RefSeq=${params.refseq}
				export ConversionFile=${params.conversionFile}
				calc_maximal_isoform.bash
				"""
		}

		filteredIsoforms
				.into() { filteredIsoformsForSingleRef; }

		// Use the first sample in the list (sorted by name) as the reference.
		// There isn't much consistency in isoforms when comparing across samples
		filteredIsoformsForSingleRef
				.toSortedList()
				.map() { it -> it[0] }
				.into() { singleRef }
} else {
		singleRef = Channel.fromPath(params.refseq)
				.map() { it -> tuple("RefSeq", it)}
}

// Generate a SAF using filtered isoforms for nascent data
process genIsoformSAF {
		cpus 1
		memory '1GB'
		time '10m'
		tag "$prefix"
		publishDir "${params.outdir}/isoform/", mode: 'copy', pattern: "isoform_max.gtf", overwrite: true
		input:
		set val(prefix), file(isoform_max) from singleRef

		output:
		file("isoform_max.saf") into singleRefSAF

		module 'python/3.6.3'
		module 'bedtools'
		module 'subread'
		script:
		"""
		# This is where we do the filtering
		awk -v OFS='\t' '{if (\$3-\$2 > ${params.threshold}) {print \$4, \$1, \$2, \$3, \$6}}' "${isoform_max}" > isoform_max.saf
		"""
}

// We want to generate counts separately for each bam file
allBam = bamInit
		.toSortedList()

// Step 2.1 -- Generate counts independently for speed
process perSampleCounts {
		cpus 8
		memory '2 GB'
		time '30m'
		tag "$prefix"
  publishDir "${params.outdir}/counts/individual/", mode: 'copy', pattern: "counts*", overwrite: true
	input:
		each file(bam) from bamForCounts
		file(single_saf) from singleRefSAF

	output:
		set file("*_without_header"), file(bam) into countsTableIndividualForMerge
		set file("*_without_header"), file(bam) into countsTableIndividualForSplit

	module 'python/3.6.3'
	module 'bedtools'
	module 'subread'
	script:
		"""
		# Set up variables we need
		refFlag='SAF'
		refFile=${single_saf}
		export refFlag
		export refFile
		# Export array as a function to get it into the script.
		function exportArray {
	  BamFiles=("${bam}")
}
		export -f exportArray
		bash -c \"exportArray; . counts_for_deseq.bash ${bam} 1\"
		"""
}

// Separate into different counts tables to use different length
// counting methods. We want to count the full gene for nascent and
// only exons for RNA
allCounts = countsTableIndividualForMerge
		.map{a, b -> a}
		.toSortedList()

// Step 2.2
process mergeIndividualCounts {
	cpus 1
	memory '2 GB'
	time '10m'
  tag "$prefix"
  publishDir "${params.outdir}/counts/", mode: 'copy', pattern: "counts_merged.txt*", overwrite: true
	input:
		file(count) from allCounts

		output:
		file("counts_merged.txt") into countsTableForVSI

		script:
		"""
		files=(${count})
		num_files="\${#files[@]}"
		echo "\$num_files"
		merge_counts_with_join.r -i \${files[@]} -o counts_merged.txt
		"""
}

process runVSI {
		cpus 4
		memory '64 GB'
 		time '720m'
		tag "$prefix"
		publishDir "${params.outdir}/vsi/", mode: 'copy', pattern: "*", overwrite: true
		input:
		file(counts_merged) from countsTableForVSI

		script:
		"""
		python3 /scratch/Users/zama8258/virtual_spike_in/virtual_spike_in/main.py \
				--label "virtual" \
				--outdir "${params.outdir}" \
				--count-data "${counts_merged}" \
				--generate-plots
		"""
}
