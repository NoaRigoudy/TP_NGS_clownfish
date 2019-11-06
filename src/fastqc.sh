#!/bin/bash

cd /home/rstudio/disk/data
mkdir -p fastqc_data
# Creation d'un dossier pour stocker nos donnees de fastqc


list_fqc=`ls sra_data/*.fastq`

# Creation d'une liste qui contient tous les noms de fichiers .fastq 
#que l'on veut analyser



for i in $list_fqc
do
fastqc $i -o /home/rstudio/disk/data/fastqc_data
done

# Verification de qualite des fichiers fastq grace a fastqc
#en utilisant une boucle for. On stocke dans un nouveau dossier fastqc_data
