###
### use R CMD BATCH -_AT_QT.paired_matched.fastq.gz -001_AT_QT.paired_matched.fastq.gz -FILE /groups/DAC/useful_R/illumina_HiSeq_QC_args.R
###
### find ./ -iname '*001.fastq.gz' -exec R CMD BATCH -fastq.gz -001.fastq.gz -'{}' /groups/DAC/useful_R/illumina_HiSeq_QC_args.R \;
### find ./ -iname '*001.fastq.gz' -exec R CMD BATCH -fastq.gz -001.fastq.gz -'{}' /data2/pipeline_in/Runs/useful_R/illumina_HiSeq_QC_args.R \;
###
### ARG[1]=file extension to look for in all subfiles, ARGV[2]=file extenstion to pull off to get all files via a pattern search, ARGV[3]=first file in set
require(ShortRead)
illumqc<-function(filemask){
	i=1
                f = FastqStreamer(files[i], n=10000)
                fq<-yield(f)
                numcycles<-max(width(fq))
                qual<-quality(fq)
                seq<-sread(fq)
                num_reads<-length(fq)
                widfq<-width(fq)
                qual_summary<-alphabetByCycle(qual)
                seq_summary<-alphabetByCycle(seq)
                sum_len<-as.numeric(sum(widfq,na.rm=TRUE))
                while (length(fq<-yield(f))) {
                        numcycles<-max(c(numcycles,max(width(fq))))
                        qual<-quality(fq)
                        seq<-sread(fq)
                        num_reads<-num_reads+length(fq)
                        widfq<-width(fq)
                        sum_len<-sum_len+sum(widfq,na.rm=TRUE)
                        qual_summary<-qual_summary+alphabetByCycle(qual)
                        seq_summary<-seq_summary+alphabetByCycle(seq)
                }
                close(f)
                rm(fq)
	i=2
	while(i<(length(files)+1)){
		f = FastqStreamer(files[i], n=10000)
		fq<-yield(f)
		len_fq<-length(fq)
		if(len_fq>0){
			numcycles<-max(width(fq))
			qual<-quality(fq)
			seq<-sread(fq)
			num_reads<-num_reads+length(fq)
			widfq<-width(fq)
			qual_summary<-qual_summary+alphabetByCycle(qual)
			seq_summary<-seq_summary+alphabetByCycle(seq)
			sum_len<-sum_len+sum(widfq,na.rm=TRUE)
			fq<-yield(f)
			len_fq<-length(fq)
			while (len_fq>0) {
				print(length(fq))
				numcycles<-max(c(numcycles,max(width(fq))))
				qual<-quality(fq)
				seq<-sread(fq)
				num_reads<-num_reads+length(fq)
				widfq<-width(fq)
				sum_len<-sum_len+sum(widfq,na.rm=TRUE)
				qual_summary<-qual_summary+alphabetByCycle(qual)
				seq_summary<-seq_summary+alphabetByCycle(seq)
				fq<-yield(f)
				len_fq<-length(fq)
			}
		}
		close(f)
		rm(fq)
		i=i+1
	}

	avg_len<-as.integer((sum_len)/num_reads)
	mat<-qual_summary
	vals<-rownames(mat)
	vals<-paste(vals, collapse="")
	vals<-as.numeric(charToRaw(vals))
	vals<-vals-33
	mat2<-ifelse(mat == 0, NA, mat)
	mat2<-mat2*vals
	means<-(apply(mat2, 2, sum, na.rm=TRUE))/(apply(mat, 2, sum, na.rm=TRUE))
	sds<-(apply(mat2, 2, sd, na.rm=TRUE))/(apply(mat, 2, sum, na.rm=TRUE))
	highs<-apply((vals*mat2/mat2),2,max,na.rm=TRUE)
	lows<-apply((vals*mat2/mat2),2,min,na.rm=TRUE)
	bmp<-paste(filepattern, "_plot.png", sep="")
	bitmap(bmp, res=300)
	
	subtitle<-paste("number of reads: ", num_reads, " avg read length: ", avg_len, sep="") 
	par(oma=c(3,3,3,3)) 
	plot(means, ylab="Mean quality score", xlab="Cycle Number", ylim=c(0, 42))
	title(filepattern, line=3)
	title(subtitle, line=2, cex.main=0.8)
	arrows(seq(1,numcycles), means+sds, seq(1, numcycles), means-sds, length=0.01, angle=90, code=3)
	abline(h=30, col="red", lty=2)
	#text(2, 31, labels="High quality",  pos=4)

	#Nucleotide frequency profile
	seqmat<-seq_summary[c(1:4,15),]
	seqmat3<-sweep(seqmat,2,colSums(seqmat),`/`)
	scaleSeqmat3<-max(seqmat3)
	seqmat3<-seqmat3*30/max(seqmat3)

	tick1<-(as.integer(100*scaleSeqmat3/3))/100
	tick2<-(as.integer(100*scaleSeqmat3*2/3))/100
	tick3<-(as.integer(100*scaleSeqmat3))/100

	lines(1:numcycles, seqmat3[1,], type="l", col="green")
	lines(1:numcycles, seqmat3[2,], type="l", col="blue")
	lines(1:numcycles, seqmat3[3,], type="l", col="black")
	lines(1:numcycles, seqmat3[4,], type="l", col="yellow")
	lines(1:numcycles, seqmat3[5,], type="l", col="red")
	axis(4, at=c(0, 10, 20, 30), labels=c('0', tick1, tick2, tick3))
	mtext(side=4,text="Fraction of base at position",line=3)
	legend(numcycles*1.05,42.5, c("A","C","G","T","N"), cex=0.7, col=c("green","blue","black","yellow","red"),pch=16,xpd=TRUE)
	dev.off()
	
	final_data<-rbind(means,sds)
	final_data<-rbind(final_data,seqmat3)
	#final_data_col1<-cbind(name,name,name,name,name,name,name)
	#final_data<-cbind(final_data_col1,final_data)
	#final_data[,1]<-c(name,name,name,name,name,name,name)
	rownames(final_data)<-c("means","sds","A","C","G","T","N")
	csv_name<-paste(files[1], ".dtmtx.csv", sep="")
	write.table(final_data,file=csv_name,sep=",",row.names=TRUE)
	count_data<-rbind(c("read_count","average_len"),c(num_reads,avg_len))
	#colnames(count_data)<-c("read_count","average_len")
	csv_name<-paste(files[1], ".stats.csv", sep="")
	write.table(count_data,file=csv_name,sep=",",row.names=FALSE,col.names=FALSE)
}

args <- commandArgs(trailingOnly = F)
args
myargument <- args[(length(args)-2):length(args)]
myargument<-sub("-","",myargument)
myargument
allfiles<-myargument[1]
firstfile<-myargument[3]
fileEXT<-myargument[2]
filepattern<-sub(fileEXT,"",firstfile)
filepattern<-sub("^./","",filepattern)
files<-dir(pattern=filepattern)
files
fileEXT2<-paste(allfiles,"$",sep="")
files<-grep(fileEXT2,files,perl=TRUE,value=TRUE)
files

illumqc(files)

q(save="yes")
