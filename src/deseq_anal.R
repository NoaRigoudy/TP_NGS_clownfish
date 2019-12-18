#Chargement des librairies

library("tximport")
library("readr")
library("ggplot2")

dir<-"/home/rstudio/disk/data/count_data"

#Table de correspondance entre individu, SRR, type de tissu (orange,blanc)

samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)

#Creation d'une table de correspondance gene transcrit

tx2gene=read.table(file="/home/rstudio/disk/data/trinity_data/Trinity.fasta.gene_trans_map")
#on utilise la table du output de Trinity gene_trans_map


tx2gene=data.frame(tx2gene$V2, tx2gene$V1)
names(tx2gene)=c("Transcrit","Gene")
#on echange les colonnes car tximport necessite d'avoir les transcrits en premiere colonne

files <- file.path(dir,paste(samples$SRR,".fastq", sep=""),"quant.sf")
names(files) <- samples$SRR
# creation d'une variable condition qui contient:
#dir: le chemin de travail
#samples$SRR (il faut rajouter .fastq et enlever les separateurs car les noms de dossiers 
#du count_data contiennent .fastq)
#quant.sf: le fichier output de salmon

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
#importe les donnees de quantification dans R

library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ Indiv + Tissue)

# Creation d'une table des comptes pour l'analyse DESeq (DESeqDataSet)
#a partir des donnees de quantification.
# On utilise un design Indiv + Tissue pour prendre en compte dans notre analyse
#l'effet tissue en sachant l'individu.

# Pre-filtering des donnees 

keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]
# Permet d'eliminer les genes qui n'ont pas assez de reads, leur analyse utiliserait
#de la memoire pour rien car leur analyse ne serait pas possible

# Assigner des niveaux de facteurs

dds$Tissue <- factor(dds$Tissue, levels = c("O","W"))
#on precise que orange (O) est la valeur "temoin"

#DESeq Analysis

dds <- DESeq(dds)
res <- results(dds)
summary(res)

res2<- results(dds, alpha=0.05)
summary(res2)
head(res2)
# summary des resultats de l'analyse DESeq seulement pour les points avec des p-adj<0.05 
# Analyse DESeq de notre table des comptes: 
#baseMean: valeurs d'expression normalisees
#pvalue et padj: significativite de la difference d'expression de chaque gene entre les 2 tissus

#Log fold change shrinkage
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Tissue_W_vs_O", type="apeglm")
summary(resLFC)

sum(resLFC$padj < 0.05, na.rm=TRUE)
#somme des genes qui ont une padj<0.05, on trouve 501

geneID_Trinity = read.table("/home/rstudio/disk/data/count_data/Aoce_transdecoder.stegastes.tsv", h=T)
geneID_Trinity$id=gsub("_i.*","",geneID_Trinity$Aoce)

resLFC$id=rownames(resLFC)
table_resLFC=data.frame(resLFC)

Goku = merge(Table_resLFC, geneID_Trinity,by.x="id",by.y="id" )
head(Goku)

#Search for Saiyan

Goku[Goku$Steg_ref=="ENSSPAG00000013419|si:ch211-256m1.8", ]

#Count data for Saiyan

txi$counts["TRINITY_DN20767_c0_g1",]

#plot les counts
saiyan_plot <- plotCounts(dds, gene="TRINITY_DN20767_c0_g1", intgroup="Tissue", 
                          returnData=TRUE)
library("ggplot2")
ggplot(saiyan_plot, aes(x=Tissue, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#Count data for TOP 10 genes compared to paper
#we get the gene names on ensembl and thanks to annote_fasta_from_blast.R script

#gene fhl2a
blast_stegastes<-read.table("/home/rstudio/disk/data/blast_data/blast_stegastes", h=F) 
blast_stegastes[blast_stegastes$V2 %in% "ENSSPAG00000006760.1fhl2b" ,]
#on trouve l'id Ensembl en cherchant les orthologues du gene du zebrafish chez stegastes 
#fhl2a dans le blast: l'id Trinity est 

# PC Analysis

#Transformation des donnees de compte pour la PCA

rld <- rlog(dds, blind=FALSE)
#ici on utilise rlog qui permet de prendre considere des differences a priori
#dans nos comptes selon le tissu
head(assay(rld), 3)

#PCA data and plot

PCA_data1<-plotPCA(rld, intgroup=c("Tissue", "Indiv"), returnData=T)
PCA_data2<-plotPCA(rld, intgroup=c("Tissue"))
PCA_data3<-plotPCA(rld, intgroup=c("Indiv"))

percentVar <- round(100 * attr(PCA_data1, "percentVar"))
#recupere les pourcentages de variation expliques par PC1 et PC2

ggplot(PCA_data1, aes(PC1, PC2, color=Tissue, shape=as.character(Indiv))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#Permet de mettre en forme le graphe de PCA en differenciant les tissus par couleur
#et les indivs par forme
#Attention, Indiv a transformer par as.character car ils sont numerotes et reconnus 
#comme variable continue

#MA Plot

plotMA(res, ylim=c(-5,5))
plotMA(resLFC, ylim=c(-5,5))
#Distribution des variances 














