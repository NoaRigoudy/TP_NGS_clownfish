#!/bin/bash

# Create a working directory:
data="/home/rstudio/disk/data"
mkdir -p $data    
cd $data

# -p permet de créer un dossier "data"" s'il n'existe pas deja 


# Create a directory where the data will be downloaded
mkdir -p sra_data
cd sra_data


# Make a list of SRR accessions:
SRR="SRR7591064
SRR7591067
SRR7591065
SRR7591066	
SRR7591068
SRR7591069
"

for i in $SRR 
do
fastq-dump $i

# Ne pas oublier le $ pour pouvoir rentrer dans SRR (i dans boucle)
# fastq-dump permet de récupérer les fichiers fastq depuis donnees de l'article 
# -X 4: permet de selectionner 4 spots = coordonnées X,Y d'un read sur la plaque 
#Illumina dans nos fichiers SRR




awk  '{ if (NR%2 == 1 ) {gsub("\\.","_");print $1"/1"}  else  { print $0}}' $i".fastq" > $i"_2.fastq" 



# awk permet de renommer noms de sequences en remplacant "." par "_" et ajouter /1 a la fin
# nom de sequence doit finir par /1 car trinity peut fonctionner pour des sequencages en 
#single read (/1) ou double read (deux sequences par read: /1 et /2)



rm $i".fastq"
mv $i"_2.fastq" $i."fastq"

# remplace les fichiers .fastq de depart par les fichiers de sequence renommes ("_2.fastq"")



done


