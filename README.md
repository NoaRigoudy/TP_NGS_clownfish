# NGS_Practicals_Clownfish

## Table of Contents

* Introduction
* Installation
* Usage 
  * Step 1: Data upload 
  * Step 2: Data quality control
  * Step 3: Data assembly 
  * Step 4: Data annotation 
  * Step 5: Transcript expression quantification
  * Step 6: Differential expression (DE) analysis
  * Step 6bis: Visualization of our data
  * Step 7: Searching for Saiyan, our gene of interest
* License
* Feedback/Issues
* Citation


## Introduction

This project reproduces the analysis of the clownfish transcriptome published by *Salis et al. 2019* and tries to identify *saiyan*, a candidate gene for iridophore expression in white stripes, using RNA Seq data.
[link to article] (https://onlinelibrary.wiley.com/doi/full/10.1111/pcmr.12766)

The analysis of the gene expression data uses 6 of the original datasets. Samples are taken from 2 different tissue-types (white or orange) in 3 individuals (3x2).
We analyse:

* 3 "white tissue" samples: SRR7591064, SRR7591067, SRR7591065.
* 3 "orange tissue" samples: SRR7591066, SRR7591068, SRR7591069.


## Installation

NGS_Practicals_Clownfish has the following dependencies:

Required dependencies:

* FastQC
* BWA-MEM
* SAMtools
* Integrative Genomics Viewer (IGV)
* Genome Analysis ToolKit (GATK)
* Picard Suite
* R
* Bioconductor
* DESeq


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

We want to identify the genes present in our data using Blastn. We choose to work on nucleotide sequences we can distinguish finer differences in expression; proteins being too conserved. 

Run **blast.sh** script . 

The script identifies homologies between the ORFs found in our assembled transcripts and the known transcripts of S.partitus. We use **nBlast** on the .cds file produced by **Transdecoder.LongOrfs** in "1." to blast it against the reference transcriptome downloaded in "2.". 

There are 3 steps to this script:

* build a Blastn reference database with the **makeblastdb** command, using the Stegastes partitus genome previously downloaded.
* blast the .fasta against the Blastn reference database to find homologies using the **blastn** command. 
* create a .csv annotation table running the **annote_fasta_from_blast.R** script. 
### Step 5: Transcript expression quantification

We choose to use the **Salmon** tool to quantify the expression levels of our transcripts as it is fast, accurate and bias-aware. 

Run **salmon.sh** script.

This script also has 2 parts and uses the **Salmon** program.

#### 1. Create an index transcriptome

First, we create an index or target transcriptome that Salmon uses to quasi-map RNA-Seq reads during the quantification. This uses the **salmon index ** command which needs the .fasta file as input and outputs a index folder, used by the 2nd part of Salmon.

#### 2. Quantify transcript expression in assembled data

Then, we quantify transcript-level abundances using the **salmon quant** command. It needs for input the raw data SRRxxx.fastq files and the transcript index built in "1.". It outputs a count data folder for each SRRxxx.fastq file, each containing a quant.sf file. 

These quant.sf files are then used to establish count tables. 


### Step 6: Differential expression (DE) analysis

Run **deseq_anal.sh** script. The **DESeq2** R package is needed. 

#### 1. Preparation of our data for analysis

For our DE analysis, and to use the DESeq function, three data sets are needed in table format: 

* a gene-transcript correspondance table, built using the *Trinity.fasta.gene_trans_map* output from our Trinity analysis
* the quantification output from Salmon, imported using the **tximport** command
* a metadata table as a "samples.txt" file, containing the following columns: individual_ID, SRR, tissue type (orange,white)

A count table is built using the **DESeqDataSetFromTximport** command. We here choose to include Individual + Tissue in our design to take into account the effect of each tissue (orange, white) given the individual in our analysis. 

#### 2. Pre-filtering of data and reference level specification

We pre-filter our data by eliminating all genes that have less than 10 counts from our Txi table. This saves up on memory use as performing a DE analysis for so little counts is useless.

The reference level for the Tissue parameter as Orange is set using the **factor** command. 

#### 3. DESeq analysis and Log-Fold Change shrinkage

DE analysis is performed using the **DESeq** command. 

Results can be retrieved using the **results** command, and visualized using the **summary** command. 

We perform a Log-Fold Change shrinkage on our DESeq data using the **lfcShrink** command. 
This is essential as it lowers the noise generated by genes with a low count value. Indeed, extreme differences can appear between tissues only as an artefact of low counts and may bias the analysis 
*ex: 10-->20, interpreted as x2 fold change, but in reality only the by-product of stochasticity*


### Step 6bis: Visualization of our data

We want to visualize our data to have a general overview and better apprehend our analyses. 

Run **deseq_anal.sh** script.

#### 1. PCA Analysis

A Principal Component Analysis is performed using the **plotPCA** command. 

Data must first be log-transformed using the **rlog** command. We choose to use **rlog** as it allows to consider *a priori* tissue-specific differences in gene expression levels.

#### 2. MA Plot

LFC against counts are plotted using the **plotMA** command. This enables to visualize genes that are significantly differentially expressed (p<0.05).


Both MA and PC analyses can be plotted using the **ggplot** command. 


### Step 7: Searching for Saiyan, the gene of interest

Run **deseq_anal.sh** script.

Creation of a master table called *Goku* containing for each gene:
* the Ensembl gene identifier for Stegastes reference
* the Trinity gene identifier
* the DESeq results 
* the blastn analysis data (e-value and bitscore).

Saiyan gene expression levels are retrieved using its Ensembl identifier in the Goku table. Previous search for the Ensembl ortholog is performed on the Ensembl online database. 

Count values are retrieved and then plotted using the **plotCounts** command. This enables to visualize the expression levels in each tissue and qualitatively assess differences. 
Significance of these differences are assessed thanks to the DESeq results. 

This can be done for all top 10 genes referenced in the original paper. 
*ex: gene fhl2a*

##License

Bio-RNASeq is free software, licensed under GPLv3.

## Feedback/Issues
Please report any issues to the issues page or email *noa.rigoudy@ens-lyon.fr*. If the issue calls for specific expertise, please contact the masterminds behind the project: *corentin.dechaud@ens-lyon.fr* and *marie.semon@ens-lyon.fr*.


## Citation

*P. Salis, T. Lorin, V. Lewis, C. Rey, A. Marcionetti, M.L. Escande, N. Roux, L. Besseau, N. Salamin, M. Semon, et al., Developmental and comparative transcriptomic identification of iridophore contribution to white barring in clownfish,Pigment Cell Melanoma Res, 32 (2019), pp. 391-402*









