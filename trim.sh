#!/bin/bash

#SBATCH -c 4                # number of core to be used
#SBATCH -t 0-06:00          # estimated run-time in D-HH:MM
#SBATCH -p short            # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=10000         # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB

# Get sample name
sample=${PWD##*/}

# get paths, tools

eautils=/path/fastq_tool/ExpressionAnalysis-ea-utils-bd148d4/clipper/ # installed on 21.8.20
$eautils/fastq-mcf -h
fastqc=/path/to/fastqc/FastQC_v0.11.9/FastQC/fastqc
$fastqc -version


# quality control via fastQC, plots for untrimmed sequences
date
echo "Quality control for untrimmed sequences of $sample..."
$fastqc -t 2 *.f*q.gz
echo "Done."
# move to folder
mv *html ../../fastqc/
rm *fastqc.zip

# trim sequences
date
echo "Trim sequences for $sample..."
$eautils/fastq-mcf ../../illu_ad_sel.fa *1.f*q.gz *2.f*q.gz -o ../../trim/$sample.R1.fq -o ../../trim/$sample.R2.fq -q 20 -H 
# zip trimmed sequences
cd ../../trim
pigz -p4 $sample.R1.fq
pigz -p4 $sample.R2.fq
echo "Done."

# quality control via fastQC, plots
date
echo "Quality control for trimmed sequences of $sample..."
$fastqc -t 2 $sample.R1.fq.gz $sample.R2.fq.gz
echo "Done."
date

# move to folder
mv $sample*html ../fastqc/
rm $sample*fastqc.zip




