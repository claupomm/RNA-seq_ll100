#!/bin/bash

#SBATCH -c 14               # number of core to be used
#SBATCH -t 0-01:00          # estimated run-time in D-HH:MM
#SBATCH -p ultrashort       # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=100000        # Memory pool for all cores (see also --mem-per-cpu); mem=10000 # memory 10GB


# Get sample name
sample=${PWD##*/}

# get paths, tools
star=/path/to/star/STAR-2.7.5c/bin/Linux_x86_64/STAR
genome=/path/to/star/index/Genomes/gencode/35/hs/star_chr # genome with chr and scaffolds but without haplotypes/alleles as recommended by STAR manual 2.7

date
echo "Alignment without joining sequences..."
$star --outFileNamePrefix ../../star/$sample. --genomeDir $genome --runThreadN 14 --readFilesIn ../../trim/$sample.R*.gz --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate --outSAMmode Full --outSAMmapqUnique 60 --outFilterMismatchNoverLmax 0.1 --readFilesCommand zcat  --limitBAMsortRAM 35000000000 > ../../bam/$sample.s.bam
echo "Done."
date

# remove unneeded folder
rm -r ../../star/$sample._STARtmp
