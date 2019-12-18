library("DESeq2")

dir<-"/home/rstudio/disk/data/count_data"

#Import data
countTable = read.table("/home/rstudio/disk/data/count_data/table_clown.csv", h=T)
geneID_Trinity = read.table("/home/rstudio/disk/data/count_data/Aoce_transdecoder.stegastes.tsv", h=T)
samples = read.table("/home/rstudio/disk/data/count_data/samples.txt", h=T)

#get rid of the end of the names of gene in geneID_Trinity
#Indeed, annoying to have names isoforme et p
geneID_Trinity$id=gsub("_i.*","",geneID_Trinity$Aoce)

#Make table where genes in countTable get names of stegastes' genes
DarkSocks = merge(countTable, geneID_Trinity,by.x=0,by.y="id" )

#Code to find why more gene from blast than in Darksocks
#geneID_Trinity$id[!geneID_Trinity$id%in%row.names(countTable)]
#"TRINITY_DN7296_c1_g1"
#ENSSPAG00000022656 TRINITY_DN7296_c1_g1_i1.p1 1.35e-73      278 TRINITY_DN7296_c1_g1

#Keep only info on counts and gene names frm blast
DarkSocksClean=DarkSocks[,2:7]
DarkSocksClean=round(DarkSocksClean)
row.names(DarkSocksClean)=DarkSocks$Steg_ref


#Make dataset for DESeq2 from out coutTable (and not from output trinity, but the same)
dds <- DESeqDataSetFromMatrix(countData = DarkSocksClean,
                              colData = samples,
                              design = ~ Indiv + Tissue)
dds<-DESeq(dds)
res<-results(dds)
head(res)

#Log fold change shrinkage
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Tissue_W_vs_O", type="apeglm")
head(resLFC)

#Ensembl orthologue saiyan chez Stegastes 
#ENSSPAG00000013419

# Search for Saiyan

#Count data for Saiyan

DarkSocksClean["ENSSPAG00000013419", ]
#On trouve saiyan dans notre table des comptes

saiyan_plot <- plotCounts(dds, gene="ENSSPAG00000013419|si:ch211-256m1.8", intgroup="Tissue", 
                          returnData=TRUE)
library("ggplot2")
ggplot(saiyan_plot, aes(x=Tissue, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


#DE for Saiyan

resLFC["ENSSPAG00000013419", ]


##A FAIRE

#chercher les genes differentiellement exprimes dans clownfish, creer liste de ces genes, 
#gorilla sur ces genes pour trouver les fonctions. Utiliser noms Ensembl des genes