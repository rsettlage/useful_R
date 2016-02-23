##summarizes all *.stats.csv files in a directory
##R CMD BATCH -file_output_prefix -sample_plus_barcode_cols -lane_cols -read_cols /groups/DAC/useful_R/lane_summary.R
##sample plus barcode cols is a period separated list, ie 1.2.3
require(beeswarm)

args <- commandArgs(trailingOnly = F)
args
myargument <- args
myargument<-sub("-","",myargument)
myargument
file_prefix<-myargument[7]
sample_cols<-myargument[8]
lane_col<-myargument[9]
read_col<-myargument[10]

file_prefix
lane_col

##get files to summarize and split the file names to get at lane assignments
files<-dir(pattern="stats.csv")
ff_list<-strsplit(files,"_")
df_file_fields<-data.frame(matrix(unlist(ff_list,use.names=FALSE),nrow=length(ff_list),byrow=T),stringsAsFactors=FALSE)

sample_name_col<-as.numeric(unlist(strsplit(sample_cols, ".", fixed = TRUE)))
df_file_fields2<-do.call(paste,c(df_file_fields[,sample_name_col],sep="_"))
df_file_fields2<-cbind(df_file_fields2,df_file_fields[,as.numeric(lane_col)])
df_file_fields2<-cbind(df_file_fields2,df_file_fields[,as.numeric(read_col)])

colnames(df_file_fields2)<-c("sample","Lane","Read")
df_files<-cbind(files,df_file_fields2)
df_files<-data.frame(df_files)

##now get all the data from the .csv files, remember they have a header
temp<-c("1","1")
x<-as.character(df_files$files)
for(i in 1:length(x)){
	temp<-rbind(temp,read.csv(x[i],header=TRUE))
}
temp<-temp[2:nrow(temp),]
df_files2<-cbind(df_files,temp)
df_files2$read_count<-as.numeric(df_files2$read_count)

##get row read sums
read_sums<-c("1","1")
z<-data.frame(table(factor(df_files2$Lane)))
zz<-as.character(z[,1])
for(i in 1:length(zz)){
	zz_lane<-zz[i]
	zz_sums<-sum(df_files2$read_count[df_files2$Lane==zz[i]])/(nrow(data.frame(table(factor(df_files2$Read)))))
	read_sums<-rbind(read_sums, c(zz_lane,zz_sums))
}
read_sums<-read_sums[2:nrow(read_sums),]
read_sums[2]<-as.integer(as.numeric(read_sums[2])/1000000)
read_sums[2]<-paste(read_sums[2],"M",sep="")
lane_sums<-paste(read_sums[1],read_sums[2],sep="=")

png_name=paste(file_prefix,"counts_by_sample_by_read.png",sep="_")
png(png_name)
boxplot(read_count~Read,data=df_files2,cex=0.5, main="Read counts by sample and read, colored by sample", ylab="Counts")
beeswarm(read_count~Read,data=df_files2,pwcol=df_files$sample, pch=18, cex=1.3, add=TRUE)
legend('topleft', lane_sums, cex=0.7)
mtext(side=1,text=file_prefix,line=3,cex=1)
dev.off()

summary_table<-df_files2
summary_table<-summary_table[summary_table$Read=="R1",]
csv_name=paste(file_prefix,"counts_by_sample_by_read.csv",sep="_")
write.table(summary_table, file=csv_name, sep=",",row.names=F)
