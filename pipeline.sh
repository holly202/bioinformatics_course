#!/bin/bash

#1 trimmomatic
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  $1 $2 \
  -baseout ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50

#2 fastqc on trimmed data
fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P \
	/home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P

mkdir ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads

mv ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/*fastqc* ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads/

#3 alignment

#3.1 index the reference genome
bwa index ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz

#3.2 alignment with BWA on trimmed data
mkdir ~/ngs_course/dnaseq_pipeline/data/aligned_data

bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50 ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P > ~/ngs_course/dnaseq_pipeline/data/aligned_data/aligned_data.sam

#3.3 convert to bam, sort and index
samtools view -h -b ~/ngs_course/dnaseq_pipeline/data/aligned_data/aligned_data.sam > ~/ngs_course/dnaseq_pipeline/data/aligned_data/aligned_data.bam

samtools sort ~/ngs_course/dnaseq_pipeline/data/aligned_data/aligned_data.bam > ~/ngs_course/dnaseq_pipeline/data/aligned_data/aligned_data_sorted.bam

samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/aligned_data_sorted.bam

#4 duplicate marking
picard MarkDuplicates I=~/ngs_course/dnaseq_pipeline/data/aligned_data/aligned_data_sorted.bam O=~/ngs_course/dnaseq_pipeline/data/aligned_data/aligned_data_sorted_marked.bam M=~/ngs_course/dnaseq_pipeline/data/aligned_data/marked_metrics.txt

samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/aligned_data_sorted_marked.bam

#5 post-alignment read filtering
samtools view -F 1796  -q 20 -o ~/ngs_course/dnaseq_pipeline/data/aligned_data/sorted_filtered.bam ~/ngs_course/dnaseq_pipeline/data/aligned_data/aligned_data_sorted_marked.bam

samtools index ~/ngs_course/dnaseq_pipeline/data/aligned_data/sorted_filtered.bam

#6 variant calling
zcat ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa

samtools faidx ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa

freebayes --bam ~/ngs_course/dnaseq_pipeline/data/aligned_data/sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq_pipeline/results/freebayes.vcf

bgzip ~/ngs_course/dnaseq_pipeline/results/freebayes.vcf

tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/freebayes.vcf.gz

#7 variant filtering
/usr/lib/vcflib/bin/vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/ngs_course/dnaseq_pipeline/results/freebayes.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/filtered.vcf

bedtools intersect -header -wa -a ~/ngs_course/dnaseq_pipeline/results/filtered.vcf -b ~/ngs_course/dnaseq_pipeline/data/chr22.genes.hg19.bed > ~/ngs_course/dnaseq_pipeline/results/bedtools_filtered.vcf

bgzip ~/ngs_course/dnaseq_pipeline/results/bedtools_filtered.vcf

tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/bedtools_filtered.vcf.gz

#8 annotation
~/annovar/convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq_pipeline/results/bedtools_filtered.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/bedtools_filtered.avinput

~/annovar/table_annovar.pl ~/ngs_course/dnaseq_pipeline/results/bedtools_filtered.avinput humandb/ -buildver hg19 -out ~/ngs_course/dnaseq_pipeline/results/bedtools_filtered -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout


