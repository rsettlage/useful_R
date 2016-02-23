###
### use R CMD BATCH -file illumina_QC_args.R
###
### find ./ -iname '*.gz' -exec R CMD BATCH -'{}' /groups/DAC/useful_R/illumina_PacBio_QC_args.R \;
###
require(ShortRead)
illumqc<-function(filename){
	print(filename)
	f = FastqStreamer(filename, n=100000)
		fq<-yield(f)
		numcycles<-max(width(fq))
		qual<-quality(fq)
		seq<-sread(fq)
		num_reads<-length(fq)
		widfq<-width(fq)
		qual_summary<-alphabetByCycle(qual)
		seq_summary<-alphabetByCycle(seq)
		sum_len<-sum(widfq,na.rm=TRUE)
	while (length(fq<-yield(f))) {
		numcycles<-max(c(numcycles,max(width(fq))))
		qual<-quality(fq)
		seq<-sread(fq)
		num_reads<-num_reads+length(fq)
		widfq<-width(fq)
		sum_len<-sum_len+sum(widfq,na.rm=TRUE)
		if(ncol(qual_summary)>max(width(fq))){
			q<-alphabetByCycle(qual)
			qq<-matrix(0,nrow=nrow(qual_summary),ncol=(ncol(qual_summary)-max(width(fq))))
			qqq<-cbind(q,qq)
			qual_summary<-qual_summary+qqq
			s<-alphabetByCycle(seq)
			ss<-matrix(0,nrow=nrow(seq_summary),ncol=(ncol(qual_summary)-max(width(fq))))
			sss<-cbind(s,ss)
			seq_summary<-seq_summary+sss
		} else if(ncol(qual_summary)<max(width(fq))){
			q<-qual_summary
			qq<-matrix(0,nrow=max(width(fq)),ncol=(max(width(fq))-ncol(qual_summary)))
			qqq<-cbind(q,qq)
			qual_summary<-qqq+alphabetByCycle(qual)
			s<-seq_summary
			ss<-matrix(0,nrow=max(width(fq)),ncol=(max(width(fq))-ncol(qual_summary)))
			sss<-cbind(s,ss)
			seq_summary<-sss+alphabetByCycle(seq)
		} else {
			qual_summary<-qual_summary+alphabetByCycle(qual)
			seq_summary<-seq_summary+alphabetByCycle(seq)
		}
	}
	close(f)
	rm(fq)
	
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
	##means<-(apply(mat2, 2, sum)/num_reads)
	##sds<-(apply(mat2, 2, sd)/num_reads)

	bmp<-paste(filename, "_plot.png", sep="")
	bitmap(bmp, res=300)
	
	subtitle<-paste("number of reads: ", num_reads, " avg read length: ", avg_len, sep="") 
	par(oma=c(3,3,3,3)) 
	ylimit<-42
	plot(means, ylab="Mean quality score", xlab="Cycle Number", ylim=c(0, ylimit))
	title(filename, line=3)
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
	csv_name<-paste(filename, ".csv", sep="")
	write.table(final_data,file=csv_name,sep=",",row.names=TRUE)
}

args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
fq_file <- sub("-","",myargument)
fq_file<-sub("./","",fq_file)

print(fq_file)
illumqc(fq_file)

q(save="no")
