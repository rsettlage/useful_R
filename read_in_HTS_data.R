datalist <- list() 
count_files <- list.files(pattern="\\.bam.sam.out$") 
temp<-NULL
tempname<-NULL

for(file in count_files) {
	stem <- gsub("\\_AT_QT.fastq.gz_THO2_SS_accepted_hits.bam.sorted.bam.sam.out$","",file)
	datalist[[stem]] <- read.delim(file,stringsAsFactors=F,header=T,sep="\t")
}
rm(stem)



temp<-read.delim(count_files[1],row.names=1,header=F,stringsAsFactors=F)
tempname<-sub("_AT_QT.fastq.gz_THO2_SS_accepted_hits.bam.sorted.bam.sam.out","",count_files[1])
for(i in 2:length(count_files)){
	temp1<-read.delim(count_files[i],row.names=1,header=F,stringsAsFactors=F)
	temp<-cbind(temp,temp1)
	tempname<-append(tempname,sub("_AT_QT.fastq.gz_THO2_SS_accepted_hits.bam.sorted.bam.sam.out","",count_files[1]))
}
colnames(temp)<-tempname
all_data<-temp
head(all_data)

rm(temp)
rm(tempname)
rm(temp1)
rm(file)
rm(i)
