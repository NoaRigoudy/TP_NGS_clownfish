#!/bin/bash

cd /home/rstudio/disk/data
# On se place dans le dossier data

multiqc fastqc_data -o /home/rstudio/disk/data/multiqc_data

# On effectue un multiqc sur l'ensemble des fichiers contenus
#dans le dossier fastqc_data et on les stocke dans un nouveau dossier multiqc_data


