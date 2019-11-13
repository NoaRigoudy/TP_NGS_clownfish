# NGS_Practicals_Clownfish

## Table of Contents

* [Introduction]
* [Set-Up] 
* [Usage] 
* [Step 1: Data upload] 
* [Step 2: Data quality control]
* [Step 3: Data assembly] 
* [Step 4: Data annotation] 
* [Step 5: Transcript expression quantification]


## Introduction

This project reproduces the analysis of the clownfish transcriptome published by *Salis et al. 2019* and tries to identify *saiyan*, a candidate gene for iridophore expression in white stripes, using RNA Seq data.
[link to article] (https://onlinelibrary.wiley.com/doi/full/10.1111/pcmr.12766)

The analysis of the gene expression data uses 6 of the original datasets. Samples are taken from 2 different tissue-types (white or orange) in 3 individuals (3x2).
We analyse:

* 3 "white tissue" samples: SRR7591064, SRR7591067, SRR7591065.
* 3 "orange tissue" samples: SRR7591066, SRR7591068, SRR7591069.




## Set-up


## Usage

### Step 1: Data upload

Dataset is retrieved from NCBI-SRA online database, in .fastq format under accession numbers SRRxxx.

Run **sra_download.sh** script. 

The script creates a working directory where the data will be downloaded and retrieves .fastq files using **fastq-dump** command. 

**awk** command enables to rename imported .fastq files so that file names end with /1 (important for later use of **Trinity** command). 

### Step 2: Data quality control

Run **fastqc.sh** script. 

The script runs quality analyses on all downloaded .fastq files using **fastqc** command. 

Run **multiqc.sh** script to compile and view all analyses in a single report. 

### Step 3: Data assembly

We use a *de novo* assembly tool as we do not have any reference genome for clownfish. We choose to use Trinity as it remains one of the most efficient tools for assembly-first techniques on RNA-Seq data. Indeed, Trinity has:
* a high sensitivity, 
* a high accuracy rate for short reads and full-length transcript reconstruction,
* and robustly identifies alternative splicing patterns and differenciate transcripts of paralogous genes. 

Run **trinity.sh** script. 

The script assembles all reads of the RNA Seq data and produces a .fasta file corresponding to the assembled data using the **Trinity** command. *see comments to understand how Trinity works*

Statistics are run to check assembled data quality using the **TrinityStats.pl** command. Output is a .txt file where contig length distribution is assessed. 

### Step 4: Data annotation

#### 1. Identify coding regions in our assembled transcripts

Run **transdec.sh** script. 

This script has 2 parts and uses the **Transdecoder** program. 

First, it finds the longest ORFs (Open Reading Frames) in our assembled transcripts, contained in the .fasta file created earlier by Trinity, using the **Transdecoder.LongOrfs** command.

*Note the .cds file in the output which is the nucleotide sequences for coding regions of the final candidate ORFs It is used for annotation in "3."*

Then, it predicts the likely coding regions using the **Transdecoder.Predict** command.

#### 2. Find and retrieve a "reference" genome

There is no reference genome for the clownfish. We therefore choose to work with the transcriptome of a closely-related fish species *Stegastes partitus* in order to annotate our assembled data. 

Run **getstegastes.sh** script. 

The script retrieves the *S.partitus* transcriptome dataset from the Biomart Ensembl dataset using the **wget** command. 

#### 3. Identify ORFs with homologies to known transcripts (S.partitus) 

Run **blast.sh** script . 

The script identifies homologies between the ORFs found in our assembled transcripts and the known transcripts of S.partitus. We use **nBlast** on the .cds file produced by **Transdecoder.LongOrfs** in "1." to blast it against the reference transcriptome downloaded in "2.". 

### Step 5: Transcript expression quantification

We choose to use the **Salmon** tool to quantify the expression levels of our transcripts as it is fast, accurate and bias-aware. 

Run **salmon.sh** script.

This script also has 2 parts and uses the **Salmon** program.

#### 1. Create an index transcriptome

First, we create an index or target transcriptome that Salmon uses to quasi-map RNA-Seq reads during the quantification. This uses the **salmon index ** command which needs the .fasta file as input and outputs a index folder, used by the 2nd part of Salmon.

#### 2. Quantify transcript expression in assembled data

Then, we quantify transcript-level abundances using the **salmon quant** command. It needs for input the raw data SRRxxx.fastq files and the transcript index built in "1.". It outputs a count data folder for each SRRxxx.fastq file, each containing a quant.sf file. 

These quant.sf files are then used to establish count tables. 

