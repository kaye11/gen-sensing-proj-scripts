
library(plyr)
library(dplyr)


multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T, sep=";")})
}

mydata=multmerge("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/raw track data/no dSi")

mergeddata <- ldply(mydata, data.frame)


#regex patterns

pattern <- "(\\w{1})(\\d{1})(-)(\\w+)(\\s)(\\w+)(-)(\\d+)"
well <- gsub(pattern,'\\1\\2',mergeddata$video)
wellvid <- gsub(pattern,'\\1\\2\\7\\8',mergeddata$video)
treatment <- gsub(pattern,'\\4\\5\\6',mergeddata$video)

mergeddata2=cbind(mergeddata, well, wellvid, treatment)

mergeddata2$ID=as.factor(paste(mergeddata2$wellvid, mergeddata2$A, sep="-"))

mergeddata2$Vlog=log(mergeddata2$V+1)

mergeddata2$time=as.factor(mergeddata2$T)

final.data <- with(mergeddata2, mergeddata2[order(ID, T),])

write.table (final.data, "d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data raw.csv", 
             sep=";", col.names=T, row.names=F)

#subset into control and phosphate

control.rawtrack <- subset (final.data, final.data$treatment=="control bead", )
Pbead.rawtrack <- subset (final.data, final.data$treatment=="P bead", )

#take ID, T, X and Y coords only per treatment

rawtrack <- final.data[,c("ID","T","X", "Y")] #whole data

control.rawtrack.simp <- control.rawtrack[,c("ID","T","X", "Y")] #per treatment control
Pbead.rawtrack.simp <- Pbead.rawtrack[,c("ID","T","X", "Y")] #per treatment Pbead


#save raw track

write.table (rawtrack, "d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/final track data raw_simplified.csv", 
             sep=";", col.names=T, row.names=F)

write.table (control.rawtrack.simp, "d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/control track data raw_simplified.csv", 
             sep=";", col.names=T, row.names=F)

write.table (Pbead.rawtrack.simp, "d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate/Pbead track data raw_simplified.csv", 
             sep=";", col.names=T, row.names=F)