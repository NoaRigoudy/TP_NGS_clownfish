cd /home/rstudio/disk/data/
mkdir -p transdec_data
cd transdec_data

# Creation d'un dossier pour stocker donnees de Transdec

wget -O stegastes.fa 'http://ensembl.org/biomart/martservice?query=
                        
<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
				
	<Dataset name = "spartitus_gene_ensembl" interface = "default" >
		<Attribute name = "coding" />
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "external_gene_name" />
	</Dataset>
</Query> '

# Telechargement de transcriptome de reference (Stegastes partitus) 
#depuis database Ensembl
