#!/bin/bash

cd /home/rstudio/disk/data
mkdir -p trinity_data


FASTQ=$(ls /home/rstudio/disk/data/sra_data/*.fastq |paste -s -d, -)

# Creation d'une liste avec tous les noms de fichers .fastq
# Utilisation de la fonction paste: -d permet de d'utiliser des delimiteurs de liste
#-s permet de coller un fichier a la fois et non en parallele
#- permet de lire l'input depuis ls




Trinity --seqType fq \
        --max_memory 50G \
        --CPU 14 \
        --single $FASTQ \
        --trimmomatic \
        --SS_lib_type R \
        --no_version_check \
        --normalize_by_read_set \
        --output /home/rstudio/disk/data/trinity_data
        
        
# Trinity permet assemblage de nos reads et prend plusieurs parametres en compte:
#seqType=type de fichier de nos reads, ici fq car ce sont des .fastq

#max_memory=quantite maximale de memoire de machine virtuelle allouee a Trinity,
#50Go choisi car 64Go disponible au total, on garde de la memoire en reserve

#CPU=nb de CPU alloues aux calculs de Trinity, 16 au total, on en utilise 14

#single=liste de sequences a assembler, ici compilees dans $FASTQ

#trimmomatic=permet de nettoyer rapidement les sequences fastq

#SS_lib_type=type de librairie pour nos reads : soit forward, soit reverse, soit unstranded

#no_version_check=empecher Trinity d'aller chercher version de Trinity la plus 
#recente sur internet

#normalize_by_read_set=permet de reduire la quantite de memoire necessaire aux calculs
#en normalisant les valeurs de coverage, reduisant les extremes = des reads a tres fort coverage
#ie. transcrits de genes fortement exprimes, qui pourraient prendre bcp de memoire pour rien. 

        
        
# Remove length and path in sequence names to avoid bug with blast ( sequence name length > 1000)
sed -re "s/(>[_a-zA-Z0-9]*)( len=[0-9]*)( path=.*)/\1/" /home/rstudio/disk/data/trinity_data/Trinity.fasta > Trinity2.fasta





