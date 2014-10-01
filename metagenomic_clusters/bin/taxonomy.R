library(data.table)
setwd(system("pwd", intern=T))

file_list <- list.files(path=".", pattern=".species.txt")

df<-NULL
for (file in file_list) {
table <- read.table(file, sep="\t", header=F, stringsAsFactors=F)
sorted<-table[order(table$V2, decreasing = T),]

taxa<-sorted$V1[1]
top<-sorted$V2[1]
total<-sum(sorted$V2)
discord<-total-top
prop<-top/total

all<-data.frame(file, taxa, top, total, discord, prop)
df<-rbind(df, all)
}

write.table(df, file="species.txt", sep="\t", col.names=T, row.names=F)


file_list <- list.files(path=".", pattern=".genus.txt")
df<-NULL
for (file in file_list) {
table <- read.table(file, sep="\t", header=F, stringsAsFactors=F)
sorted<-table[order(table$V2, decreasing = T),]

taxa<-sorted$V1[1]
top<-sorted$V2[1]
total<-sum(sorted$V2)
discord<-total-top
prop<-top/total

all<-data.frame(file, taxa, top, total, discord, prop)
df<-rbind(df, all)
}

write.table(df, file="genus.txt", sep="\t", col.names=T, row.names=F)


file_list <- list.files(path=".", pattern=".family.txt")
df<-NULL
for (file in file_list) {
table <- read.table(file, sep="\t", header=F, stringsAsFactors=F)
sorted<-table[order(table$V2, decreasing = T),]

taxa<-sorted$V1[1]
top<-sorted$V2[1]
total<-sum(sorted$V2)
discord<-total-top
prop<-top/total

all<-data.frame(file, taxa, top, total, discord, prop)
df<-rbind(df, all)
}

write.table(df, file="family.txt", sep="\t", col.names=T, row.names=F)

