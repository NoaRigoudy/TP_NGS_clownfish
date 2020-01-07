#! /bin/bash

data="/home/rstudio/disk/data"
cd $data
mkdir -p blast_data
cd $data/blast_data

cd ..
mv transdec_data/stegastes.fa blast_data
#On bouge notre fichier .fa dans le dossier "blast_data"

cd $data/blast_data



PATH="$PATH:/softwares/ncbi-blast-2.9.0+/bin"
#On ajoute "ncbi-blast-2.9.0+" dans les logiciels de la machine. Ne pas oublier /bin a la fin. 


sed -re "s/\|//g" stegastes.fa > stegastes_clean.fa
#On elimine les | des noms des sequences

#Run mart_to_blast.R pour retirer les duplicats du genome de reference



# Build reference database for blastn


makeblastdb -in /home/rstudio/disk/data/blast_data/stegastes_clean2.fa -parse_seqids -dbtype nucl -out stegastes

#-db_type: on choisit "nucl" car on blast des sequences nucleotidiques
#-parse_seqids: permet de garder les identifiants de sequence du fichier .fa originel



#Blast fasta against the ref db

blastn -db stegastes -query /home/rstudio/disk/data/trinity_data/Trinity.fasta -evalue 1e-5 -outfmt 6 -out blast_stegastes

#db: database=sequence reference contre laquelle on blast nos donnees
#query: fichier contenant sequence a blaster
#evalue: donne evalue minimum a considerer
#outfmt: format du output

#Rscript --vanilla  annote_fasta_from_blast.R $out_blast
#commande faite et exécutée par Corentin, on n'a pas eu le temps 
#de faire notre propre blastn en entier.


