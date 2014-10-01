rm(list=ls())
install.packages("NCBI2R", Sys.getenv("R_LIBS_USER"), repos = "http://watson.nci.nih.gov/cran_mirror/src/contrib" )
# install.packages("NCBI2R")
library(NCBI2R)
setwd(system("pwd", intern=T))

system(command="cat *blastxlist.txt | cut -f 2 | sort | uniq > blastxlist.txt")
system(command="cat *blastnlist.txt | cut -f 2 | sort | uniq > blastnlist.txt")

blastnlist<-read.delim(file="blastnlist.txt", sep="\t", header=F, stringsAsFactors=F)
uniquetaxa<-unique(blastnlist$V1)

taxid<-sapply(strsplit(as.character(uniquetaxa), '[|]'), "[", 1)
taxid<-as.numeric(gsub("taxid:", "", taxid))

listDF<-NULL
i=1
while (i < (length(taxid))) { 
	vector<-taxid[i:(i+199)]
	retrieved<-GetTaxInfo(vector)
	retrievedlin<-retrieved$lineage
	listDF<-rbind(listDF, retrievedlin)
	i=i+200
	}

head(listDF)

split<-by(listDF, as.factor(listDF$reqTaxId), invisible)

everythingDF<-NULL
lmdf<-NULL

for (i in 1: length(split)) {
	lmdf<-as.data.frame(split[[i]])
	kingdom<-lmdf$line.sciName[lmdf$line.rank=="superkingdom"]
	phylum<-lmdf$line.sciName[lmdf$line.rank=="phylum"]
	class<-lmdf$line.sciName[lmdf$line.rank=="class"]
	order<-lmdf$line.sciName[lmdf$line.rank=="order"]
	family<-lmdf$line.sciName[lmdf$line.rank=="family"]
	genus<-lmdf$line.sciName[lmdf$line.rank=="genus"]
	species<-lmdf$line.sciName[lmdf$line.rank=="species"]
	
	full<-lmdf$sciName[1]
	full<-gsub(" ", "_", full)
	taxon<-lmdf$taxId[1]
	designation<-paste(sep=";", kingdom,phylum,class,order,family,genus,full)
	# designation<-paste(sep=";", kingdom,phylum,class,order,family,genus)
	temp<-cbind(taxon, designation, full)

everythingDF<-rbind(everythingDF, temp)
}

uniqDF<-unique(everythingDF)

write.table(uniqDF, file="taxatable_blastn.txt", sep="\t", col.names=T, row.names=F)
system(command="perl ./bin/matchtax.pl taxatable_blastn.txt")
system(command="mv taxatable_blastn.txt.output.txt blastnlist_FINAL.txt")

