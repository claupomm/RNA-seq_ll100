# RNA-seq_ll100
The shell commands and scripts given serve for the RNA-seq analysis done for LL-100 cell lines, 100 leukemic and lymphoma cell lines + the normal cell line NC-NC (https://celldive.dsmz.de/rna/ll-100). 

Samples were prepared as follows: miRNeasy mini plus Qiagen, polyA selection, GATC library prep, mRNA-Seq, strand-specific fr-firststrand, reverse stranded, 2x151bp paired-end run, Illumina HiSeq2500.

The RNA-seq workflow includes count quantification via salmon and R analysis via DESeq2.


## Prepare the project folder and fastq files
set Project folder:
```
DIR=/path/to/Project_folder/Project_gatc_rna_ll100_s
```

create subfolder with samples:
```
mkdir $DIR/raw
mkdir $DIR/raw/sample_x $DIR/raw/sample_y
```

save or link fastq files in the corresponding sample folder:
```
cd $DIR/raw/K-562
ln -s /path/to/fastq_files/K-562_R1.fastq.gz .
ln -s /path/to/fastq_files/K-562_R2.fastq.gz .
```


## Download salmon for pseudo-alignment and counting, prepare indexing
taken from: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
```
# include salmon in your path:
export PATH=~/Programme/salmon-1.4.0_linux_x86_64/bin:$PATH

# make salmon transcript index:
cd /path/to/fasta_and_index_files/gencode/37

# download the reference transcriptome and genome for salmon index, protein coding:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.pc_transcripts.fa.gz

# whole genome:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz

# Preparing metadata: Salmon indexing requires the names of the genome targets, which is extractable by using the grep command:
grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

# Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference:
zcat gencode.v37.pc_transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa
pigz -p4 gentrome.fa

# Salmon Indexing on a slurm job submission system:
sbatch -c 16 -p short -J salmon -o salmon.out -e salmon.err <<EOF
#!/bin/sh
salmon index -t gentrome.fa.gz -d decoys.txt -p 16 -i salmon_index --gencode
EOF

# NOTE: --gencode flag is for removing extra metdata in the target header separated by | from the gencode reference. You can skip it if using other references
```


## Start pipeline, submit jobs to slurm
```
DIR=/path/to/Project_folder/Project_gatc_rna_ll100_s
cd $DIR
mkdir star bam plots fastqc tables trim counts
export PATH=~/Programme/salmon-1.4.0_linux_x86_64/bin:$PATH
SAMPLES=$(ls $DIR/raw/)

for SAMPLE in $SAMPLES
do
cd $DIR/raw/$SAMPLE
sample=${PWD##*/}
# trim sequences + quality control via fastQC, plots, 30-50min for 30M reads:
RES=$(sbatch -J $sample.1 -o $sample.1.out -e $sample.1.err ../../trim.sh)
# alignment without joining sequences, 6-30min for 30M reads:
RES2=$(sbatch --dependency=afterok:${RES##* } -J $sample.2 -o $sample.2.out -e $sample.2.err ../../align.sh)
# salmon count, >10min for 30M reads:
sbatch -c 12 -p ultrashort --dependency=afterok:${RES##* } -J $sample.3 -o $sample.3.out -e $sample.3.err <<EOF
#!/bin/sh
salmon quant -p 12 -l ISR --validateMappings -i /path/to/fasta_and_index_files/gencode/37/salmon_index -r ../../trim/> $sample.*fq.gz -o ../../counts/$sample.salmon # --gcBias recommended for PE; ISR inward stranded reverse; -fr-firststrand: -l ISR -l SR; -fr-secondstrand: -l ISF -l SF
EOF
done
```


## Checking the job status, results, remove non-neccessary files
```
cd $DIR
grep "Mapping rate" counts/*salmon/logs/*log
ls -lh $DIR/trim/*
ls -lh raw/*/*.1.err
ls -lh raw/*/*.2.out
ls -lh raw/*/*.3.out

rm raw/*/*.1.err
rm trim/*
rm raw/*/*.2.out
rm raw/*/*.3.out
```


## Info on sequencing + mapping, star statistics
```
cd $DIR
R --file="stats_seq.R"
```


## DESeq2 analysis, normalisation
Before running this R script, create a file2sample.csv data with sample assignment. Furthermore download ensembl data for annnotation.
```
R --file='analysis_deseq2.R' &> R_salmon.out &
```


