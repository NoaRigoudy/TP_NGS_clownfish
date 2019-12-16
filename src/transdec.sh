#! /bin/bash

data="/home/rstudio/disk/data"
mkdir -p transdec_data
cd $data/transdec_data



# Step 1: extract the long open reading frames

#TransDecoder.LongOrfs -t $data/trinity_data/Trinity.fasta -m 100 --gene_trans_map $data/trinity_data/Trinity.fasta.gene_trans_map -S 


# Transdecoder.LongOrfs extrait les ORFs les plus longs dans nos sequences assemblees, 
#contenues dans Trinity.fasta

#-t : input file
#-m : identifie par default des ORFs de longueur 100aa, 
#<100 donne lieu generalement a un bruit de fond trop important
#gene_trans_map: tableau qui contient les genes et leurs differents transcrits
# -S: specifie si nos transcrits sont orientes, 
#on a specifie a Trinity que notre librairie de strands etait en reverse (R) 
#donc nos transcrits sont orientes

# Output=
#transcripts.fasta.transdecoder.pep : peptide sequences for the final candidate ORFs; 
#all shorter candidates within longer ORFs were removed.

#transcripts.fasta.transdecoder.cds  : nucleotide sequences for coding regions of 
#the final candidate ORFs

#transcripts.fasta.transdecoder.gff3 : positions within the target transcripts of 
#the final selected ORFs

#transcripts.fasta.transdecoder.bed  : bed-formatted file describing O




# Step 2: Predict the likely coding regions 

TransDecoder.Predict -t $data/trinity_data/Trinity.fasta --single_best_only --cpu 10

# TransDecoder.Predict:
#-t: fichier input 
# --single_best_only: retient que le meilleur ORF par transcrit, 
#prioritisation par homologie puis longueur de l'ORF   
# -O: destination pour output file
# --cpu: nombre de CPU alloues a calculs




# (Etape suivante: Identify ORFs with homologies to known transcripts using Blastn and .cds output file)

#cf blast.sh
# On utilisera ensuite le fichier .cds pour le blaster contre des sequences nucleotidiques connues grace a Blastn



