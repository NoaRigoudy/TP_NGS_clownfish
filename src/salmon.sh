#! /bin/bash

data="/home/rstudio/disk/data"
cd $data
mkdir -p salmon_data
cd $data/salmon_data

mkdir -p transcript_index

# Index the transcriptome

#salmon index -t $data/trinity_data/Trinity.fasta -k 25 -i transcript_index

# -t: transcript fasta file
# -k: definit longueur de k-mer a utiliser pour contruire l'index,
#on considere une longueur max k=31bp et une valeur optimale de 50-66% de longueur de nos reads
#je choisis 50% ie. k=25


# Quantification

FASTQ=$(ls /home/rstudio/disk/data/sra_data/*.fastq)

# Creation d'une liste avec tous les noms de fichers .fastq

for i in $FASTQ
do

base=`basename $i`

# Creation d'une liste contenant seulement les noms des fichiers .fastq et non tout le directory
#comme ce qui est stocke dans FASTQ


salmon quant -l SR -i transcript_index -p 14 --gcBias -r $i -o "/home/rstudio/disk/data/count_data/"$base
# -l= library type, ici on est en stranded reverse (SR)
# -i= folder index
# -p= nb de coeurs a utiliser
# --gcBias= corrige biais de GC content dans nos reads
# -r= fichier fastq contenant reads
# -o: output file directory

done

#

