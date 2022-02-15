# EcoTILLING of GRANULE BOUND STARCH SYNTHASE in Arabidopsis thaliana 

## This is linked to the publication by Seung et al.: [Natural Polymorphisms in Arabidopsis Result in Wide Variation or Loss of the Amylose Component of Starch](https://pubmed.ncbi.nlm.nih.gov/31694903/)

## Supplemental Figure 6


### Obtaining sequences

Download SRR1946353 and SRR1946535 from NCBI using sra toolkit

```
~/programs/sratoolkit.2.9.6-mac64/bin/prefetch SRR1946535
~/programs/sratoolkit.2.9.6-mac64/bin/prefetch SRR1946353
~/programs/sratoolkit.2.9.6-mac64/bin/fastq-dump --split-3 --gzip SRR1946353.sra
~/programs/sratoolkit.2.9.6-mac64/bin/fastq-dump --split-3 --gzip SRR1946535.sra

```

* Gn2-3 (sample 9790) SRR1946353
* TueSB30-3 (sample 9999) SRR1946535. These data seem to be single end sequencing.

Download TAIR 10 Arabidopsis thaliana chromosomes from [arabidopsis.org](ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/) and concatenate sequences

```
cat *.fas > genome.fasta
```


### Mapping

```
bwa index genome.fasta

bwa aln genome.fasta SRR1946353_1.fastq.gz > SRR1946353_1.aln
bwa aln genome.fasta SRR1946353_2.fastq.gz > SRR1946353_2.aln
bwa aln genome.fasta SRR1946535.fastq.gz > SRR1946535.aln

bwa sampe genome.fasta SRR1946353_1.aln SRR1946353_2.aln SRR1946353_1.fastq.gz SRR1946353_2.fastq.gz >SRR1946353.sam
bwa samse genome.fasta SRR1946535.aln SRR1946535.fastq.gz > SRR1946535.sam

samtools sort -o SRR1946353.bam SRR1946353.sam
samtools sort -o SRR1946535.bam SRR1946535.sam

samtools index SRR1946353.bam
samtools index SRR1946535.bam

samtools faidx genome.fasta
samtools mpileup -BQ0 -f genome.fasta SRR1946353.bam > SRR1946353.pileup
samtools mpileup -BQ0 -f genome.fasta SRR1946535.bam > SRR1946535.pileup

```


### Compare pileup files

The java file DisplaySNPRatio.java contains the main method to read in the pileup files.

### Plot figure

run the R script plot_distance.R


