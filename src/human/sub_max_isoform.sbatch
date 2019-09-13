#!/bin/bash

echo "Submitting Multiple CalcIsoform Jobs"
dirPrefix=/scratch/Shares/dowell/ripin/PRO-seq/mapped/bams

c1="$dirPrefix"/wt_1hAs_U2OS_b1.sorted.bam
c2="$dirPrefix"/wt_1hAs_U2OS_b2.sorted.bam
c3="$dirPrefix"/wt_U2OS_b1.sorted.bam
c4="$dirPrefix"/wt_U2OS_b2.sorted.bam
p1="$dirPrefix"/dG3BP_U2OS_b1.sorted.bam
p2="$dirPrefix"/dG3BP_U2OS_b2.sorted.bam
p3="$dirPrefix"/dG3BP_1hAs_U2OS_b1.sorted.bam
p4="$dirPrefix"/dG3BP_1hAs_U2OS_b2.sorted.bam

sbatch --export	InterestFile="$c1" --job-name="C1-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$c2" --job-name="C2-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$c3" --job-name="C3-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$c4" --job-name="C4-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$p1" --job-name="P1-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$p2" --job-name="P2-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$p3" --job-name="P3-MaxIso" calc_maximal_isoform.sbatch
sbatch --export	InterestFile="$p4" --job-name="P4-MaxIso" calc_maximal_isoform.sbatch
