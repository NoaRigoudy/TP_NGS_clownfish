#Chargement des librairies

library("tximport")
library("readr")
library("ggplot2")

dir<-"/home/rstudio/disk/data/count_data"



### Step 1 : Preparation pour Analyse DESeq


## Table de correspondance entre individu, SRR, type de tissu (orange,blanc)

samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)

## Creation d'une table de correspondance gene transcrit

tx2gene=read.table(file="/home/rstudio/disk/data/trinity_data/Trinity.fasta.gene_trans_map")
#On utilise la table du output de Trinity gene_trans_map


tx2gene=data.frame(tx2gene$V2, tx2gene$V1)
names(tx2gene)=c("Transcrit","Gene")
#On echange les colonnes dans notre tableau car tximport necessite d'avoir les transcrits 
#en premiere colonne

files <- file.path(dir,paste(samples$SRR,".fastq", sep=""),"quant.sf")
names(files) <- samples$SRR
#Creation d'une variable condition qui contient:
#dir: le chemin de travail
#samples$SRR (il faut rajouter .fastq et enlever les separateurs car les noms de dossiers 
#du count_data contiennent .fastq)
#quant.sf: le fichier output de salmon


## Import des donnees de quantification dans R

txi <- tximport(files, type="salmon", tx2gene=tx2gene)



## Creation d'une table des comptes pour l'analyse DESeq (DESeqDataSet)
#a partir des donnees de quantification.

library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ Indiv + Tissue)

#On utilise un design Indiv + Tissue pour prendre en compte dans notre analyse
#de l'effet tissu en sachant l'individu.



## Pre-filtering des donnees 

keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]
#Permet d'eliminer les genes qui n'ont pas assez de reads, leur analyse utiliserait
#de la memoire pour rien car elle ne serait pas possible au vu de leur nombre faible


## Assigner des niveaux de facteurs

dds$Tissue <- factor(dds$Tissue, levels = c("O","W"))
#On precise que orange (O) est la valeur "temoin"




### Step 2 : Analyse DESeq

##Analyse 

dds <- DESeq(dds)
res <- results(dds)
summary(res)
#Summary des resultats de l'analyse DESeq

res2<- results(dds, alpha=0.05)
summary(res2)
#Summary des resultats de l'analyse DESeq seulement pour les points avec des p-adj<0.05

head(res2)
#baseMean: valeurs d'expression normalisees
#pvalue et padj: significativite de la difference d'expression de chaque gene entre les 2 tissus


## Log fold change shrinkage

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Tissue_W_vs_O", type="apeglm")
summary(resLFC)
#Summary sur les donnees shrinkees ce qui permet de reduire le bruit genere par 
#par les genes a faible nombre de copies

sum(resLFC$padj < 0.05, na.rm=TRUE)
#Somme des genes qui ont une padj<0.05, on en trouve 501



## Visualtisation des resultats


# PC Analysis

rld <- rlog(dds, blind=FALSE)
# Transformation des donnees de comptes pour la PCA
#ici on utilise rlog qui permet de prendre considere des differences a priori
#dans nos comptes selon le tissu

head(assay(rld), 3)


# PCA data and plots

PCA_data1<-plotPCA(rld, intgroup=c("Tissue", "Indiv"), returnData=T)
PCA_data2<-plotPCA(rld, intgroup=c("Tissue"))
PCA_data3<-plotPCA(rld, intgroup=c("Indiv"))
#Permet de visualiser les principaux parametres expliquant la variance observee

percentVar <- round(100 * attr(PCA_data1, "percentVar"))
#Recupere les pourcentages de variation expliques par PC1 et PC2


ggplot(PCA_data1, aes(PC1, PC2, color=Tissue, shape=as.character(Indiv))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#Permet de mettre en forme le graphe de PCA en differenciant les tissus par couleur
#et les indivs par forme
#Attention, Indiv a transformer par as.character car ils sont numerotes et reconnus 
#comme variable continue


# MA Plot

plotMA(res, ylim=c(-5,5))
plotMA(resLFC, ylim=c(-5,5))
#Permet de visualiser les genes qui sont differentiellement exprimes.




### Step 3 : Search for Saiyan


geneID_Trinity = read.table("/home/rstudio/disk/data/count_data/Aoce_transdecoder.stegastes.tsv", h=T)
geneID_Trinity$id=gsub("_i.*","",geneID_Trinity$Aoce)
#Creation d'une colonne id contenant des identifiants nettoyes pour chaque gene dans la table de correspondance entre identifiant de gene Trinity, 
#identifiant Ensembl, et e-value et bitscore du blastn. 

resLFC$id=rownames(resLFC)
table_resLFC=data.frame(resLFC)
#Creation d'une table contenant les resultats de l'analyse differentielle DeSeq


Goku = merge(Table_resLFC, geneID_Trinity,by.x="id",by.y="id" )
head(Goku)
#Creation d'une table combinant les correspondances entre identite des genes Trinity, 
#les resultats du blastn et nos valeurs pour l'analyse d'expression differentielle.
#Ceci permet d'avoir les niveaux d'expression pour chaque transcrit identifie.


Goku[Goku$Steg_ref=="ENSSPAG00000013419|si:ch211-256m1.8", ]
#On cherche les niveaux d'expression de notre gene d'interet Saiyan grace a son identifiant Ensembl 
#dans notre table Goku


## Count data for Saiyan

txi$counts["TRINITY_DN20767_c0_g1",]
#Cherche les comptes pour le gene Saiyan en particulier

#Plot des comptes

saiyan_plot <- plotCounts(dds, gene="TRINITY_DN20767_c0_g1", intgroup="Tissue", 
                          returnData=TRUE)
library("ggplot2")
ggplot(saiyan_plot, aes(x=Tissue, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


#Comparaison des donnees de comptes pour notre TOP 10 de genes differentiellement exprimes, 
# a comparer avec le papier


#On recupere les noms des genes sur Ensembl et grace au script annote_fasta_from_blast.R 

#gene fhl2a
blast_stegastes<-read.table("/home/rstudio/disk/data/blast_data/blast_stegastes", h=F) 
blast_stegastes[blast_stegastes$V2 %in% "ENSSPAG00000006760.1fhl2b" ,]
#on trouve l'id Ensembl en cherchant les orthologues du gene du zebrafish chez stegastes 
#fhl2a dans le blast: l'id Trinity est 















