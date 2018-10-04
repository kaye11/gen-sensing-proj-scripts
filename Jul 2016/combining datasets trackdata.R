
library(dplyr)
library(plyr)

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T, sep=";")})
}

mydata=multmerge("d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate automatic/binned track data")

mergeddata <- ldply(mydata, data.frame)


#regex patterns

pattern <- "(\\w{1})(\\d{1})(-)(\\w+)(\\s)(\\w+)(-)(\\d+)"
well <- gsub(pattern,'\\1\\2',mergeddata$video)
wellvid <- gsub(pattern,'\\1\\2\\7\\8',mergeddata$video)
treatment <- gsub(pattern,'\\4\\5\\6',mergeddata$video)

mergeddata2=cbind(mergeddata, well, wellvid, treatment)

mergeddata2$ID=as.factor(paste(mergeddata2$wellvid, mergeddata2$A, sep="-"))

mergeddata2$Vlog=log(mergeddata2$V+1)

mergeddata2$time=as.factor(mergeddata2$time)

mergeddata2 <- with(mergeddata2, mergeddata2[order(ID, bin, time),])

#drop outbin
final.data <- mergeddata2[mergeddata2$bin!= "outbin",]

write.table (final.data, "d:/Karen's/PhD/R program/General sensing proj/csv files/Tracking phosphate automatic/final track data.csv", 
             sep=";", col.names=T, row.names=F)