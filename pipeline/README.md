# Virtual Spike-In Pipeline Implementation

This folder contains a nextflow pipeline implementing the virtual spike-in pipeline from start to finish.
The pipeline requires the following dependencies:
- `nextflow`
- `bedtools`
- `featureCounts`
- `virtualspikein`

To use the pipeline, configure the following parameters in `config/base.conf`:
- `bamDir`: (Path) The path to a folder containing your bam files
- `refseq`: (Path) The path to a bedfile containing your reference sequence
- `conversionFile`: (Path) The path to a file containing tab delimited key/value pairs converting NCBI IDs to Gene Symbols
- `timepoint`: (Int) The timepoint of your samples in minutes
- `elongationRate`: (Int) The elongation rate of RNA polymerase in your organism in bases/min (Default: 3000 for human)
- `debugInfo`: (Bool) Whether to return additional diagnostic information for debugging
- `outdir`: (Path) The path of an existing directory to send output files to

In addition, make sure that you've read the requisite warnings in the main README about when this tool is not suitable for your data. In brief, this tool and pipeline are insufficient for your data if you are (a) running the virtual spike-in only and (b) you suspect your experimental perturbation directly targets your organism's RNA Polymerase.

This pipeline will return output in the following format:
- `{outDir}/isoform` - Reference bedfiles filtered to include only a single isoform
- `{outDir}/counts` - Counts files for each individual sample and a merged counts table for the tool to use
- `{outDir}/vsi` - Descriptive output including: size factors, uncertainty bounds, diagnostic information
