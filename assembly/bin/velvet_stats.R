rm(list=ls())
library(ggplot2)
library(reshape2)

setwd("./skinmetagenome/assembly/bin")

input.file <- "bowtiestats.txt";
output.file <- "best.txt";
table<-read.delim(file=input.file, sep="\t", header=F, stringsAsFactors=F)

label<-sapply(strsplit(as.character(table$V1), '[ :]'), "[", 1)
kmer<-as.numeric(sapply(strsplit(as.character(label), '[.]'), "[", 2))
type<-sapply(strsplit(as.character(label), '[.]'), "[", 3)
concord1<-as.numeric(sapply(strsplit(as.character(table$V1), '[ :]'), "[", 2))
discord1<-as.numeric(sapply(strsplit(as.character(table$V2), '[ :]'), "[", 2))
overall<-sapply(strsplit(as.character(table$V3), '[ :]'), "[", 2)
overall<-as.numeric(gsub("%","", overall))
totalpaired <-as.numeric(sapply(strsplit(as.character(table$V4), '[ :]'), "[", 2))
totalreads<-as.numeric(sapply(strsplit(as.character(table$V5), '[ :]'), "[", 2))

percconcord<-concord1/totalpaired
percdiscord<-discord1/totalpaired
overall<-overall/100

newtable<-data.frame(label, kmer, type, percconcord, percdiscord, overall)
attach(newtable)
names(newtable)

##############################################################
input.file2<-"all_assemstats.txt"
table2<-read.table(file=input.file2, sep="\t", header=F, stringsAsFactors=F)
names(table2)<-c("N", "sum", "max", "filename")
attach(table2)

kmer2<-as.numeric(sapply(strsplit(as.character(table2$filename), '[./]'), "[", 2))
newlabel2<-table2$filename
category<-rep("velvet", nrow(table2))

newtable2<-data.frame(table2$filename, newlabel2, category, kmer2, N, sum, max)

##############################################################
head(newtable)
head(newtable2)

bowtietable<-data.frame(newtable, totalreads)
assemtable<-newtable2
assemtable$matcher<-gsub("_metavelvet_assembly", "", assemtable$table2.filename)
assemtable$matcher<-gsub("_velvet_assembly", "", assemtable$matcher)
assemtable$matcher<-gsub(".fa", "", assemtable$matcher)
assemtable$matcher<-paste(assemtable$matcher, category, "stats.txt", sep=".")

write.table(bowtietable, file="bowtietable.txt", sep="\t", col.names=T, row.names=F)
write.table(assemtable, file="assemtable.txt", sep="\t", col.names=T, row.names=F)

both<-data.frame(assemtable[1:7], bowtietable[4:7][match(assemtable$matcher, bowtietable$label),])
write.table(both, file="both.txt", sep="\t", col.names=T, row.names=F)
both<-na.omit(both)
velvetonly<-both[both$category =="velvet",]
newlabel<-refs$Illumina.Code[match(velvetonly$V1, refs$All)]
velvetonly$V2<-sapply(strsplit(as.character(velvetonly$table2.filename), '[_]'), "[", 1)

split<-by(velvetonly, as.factor(velvetonly$V2), invisible)


#give everyone a rank
#then *2 their rank for overall alignment
#if = 0, rank as BAD. 

everythingDF<-NULL
lmdf<-NULL
for (i in 1: length(split)) {
	lmdf<-as.data.frame(split[[i]])
lmdf$Nrank<-rank(lmdf$N)
lmdf$sumrank<-rank(lmdf$sum)
lmdf$overallrank<-rank(lmdf$overall)
lmdf$score<-lmdf$Nrank + lmdf$sumrank + 2*lmdf$overallrank
lmdf<-lmdf[order(lmdf$score, decreasing=TRUE),]
	# ordered<-lmdf[order(lmdf$filename, -lmdf$score),]
	# use<-ordered[1,]
use<-lmdf[1,]
everythingDF<-rbind(everythingDF, use)
}

everythingDF$label<-paste(everythingDF$V2, "_assembly.", everythingDF$kmer2, ".fa", sep="")
write.table(everythingDF, file="best.velvet.txt", col.names=T, row.names=F, sep="\t")
write.table(as.character(everythingDF$label), file="velvet_filestoget.txt", col.names=F, row.names=F, sep="\t")

