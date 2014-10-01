install.packages("NCBI2R", Sys.getenv("R_LIBS_USER"), repos = "http://watson.nci.nih.gov/cran_mirror/src/contrib" )
# install.packages("NCBI2R")
library(NCBI2R)
setwd(system("pwd", intern=T))

blastnlist<-read.table(file="allgis.txt.output.txt", sep="\t", header=F)
taxid <-unique(blastnlist$V2)

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

uniqDF<-data.frame(unique(everythingDF))

#now, we have to match to the GIs
head(blastnlist)
head(uniqDF)
blastnlist$designation<-uniqDF$designation[match(blastnlist$V2, uniqDF$taxon)]
blastnlist$full<-uniqDF$full[match(blastnlist$V2, uniqDF$taxon)]
names(blastnlist)<-c("gi", "taxid", "designation", "full")

write.table(blastnlist, file="blastxlist_FINAL.txt", sep="\t", col.names=T, row.names=F)


