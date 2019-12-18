
install.packages("seqinr")
library(seqinr)
f=read.fasta("/home/rstudio/disk/data/blast_data/stegastes_clean.fa")

#Remove duplicates from stegastes_clean.fa (without |) 

sequences <- sapply(f, c2s)

sequences = sequences[sequences!="sequence unavailable"]
genes <- unique(names(f))

filtered_seq <- lapply(genes, function(v) {
  
  sq=sequences[names(sequences)==v]
  
  if (length(sq) != 0) {
    lsq=sapply(1:length(sq),function (x) {
      return(unname(nchar(sq[x])))
    })
    
    return(c(v,toupper(sq[which(lsq==max(lsq))][1])))
  }
})

vector = unlist(filtered_seq)
names = vector[c(TRUE,FALSE)]
sequences = as.list(vector[c(FALSE,TRUE)])


write.fasta(sequences=sequences, names=names, file.out="/home/rstudio/disk/data/blast_data/stegastes_clean2.fa")
