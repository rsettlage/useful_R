require(ShortRead)
illumqc<-function(filename){
	#sampler = FastqStreamer(filename, n=10000)
	sampler = FastqSampler(filename, n=100000)
	fq = yield(sampler)
	#fq<-readFastq(filename)
	numcycles<-max(width(fq))
	qual<-quality(fq)
	seq<-sread(fq)
	lenfq<-length(fq)
	widfq<-width(fq)
	rm(fq)
	mat<-alphabetByCycle(qual)
	vals<-rownames(mat)
	vals<-paste(vals, collapse="")
	vals<-as.numeric(charToRaw(vals))
	vals<-vals-33
	mat2<-mat*vals
	means<-(apply(mat2, 2, sum)/lenfq)
	sds<-(apply(mat2, 2, sd)/lenfq)
	name=filename
	#sections<-strsplit(filename, "\\.|\\/")[[1]]
	#for (i in 1:length(sections)) {
  	#	if (sections[i] == 'fastq') {name=sections[i-1]}
  	#} 
	bmp<-paste(name, "_plot.png", sep="")
	bitmap(bmp, res=300)

	plot(means, ylab="Mean quality score", xlab="Cycle Number", ylim=c(0, 42), main=name)
	arrows(seq(1,max(widfq)), means+sds, seq(1, max(widfq)), means-sds, length=0.01, angle=90, code=3)
	abline(h=30, col="red", lty=2)
	text(2, 31, labels="High quality",  pos=4)

	#Nucleotide frequency profile
	seqmat<-alphabetByCycle(seq)
	seqmat2<-rbind(seqmat[1:4,], seqmat[15,])
	seqmat3<-80*seqmat2/lenfq

	lines(1:numcycles,seqmat3[1,], type="l", col="green")
	lines(1:numcycles, seqmat3[2,], type="l", col="blue")
	lines(1:numcycles, seqmat3[3,], type="l", col="black")
	lines(1:numcycles, seqmat3[4,], type="l", col="yellow")
	lines(1:numcycles, seqmat3[5,], type="l", col="red")
	axis(4, at=c(0, 10, 20, 30), labels=c('0', '0.125', '0.25', '0.375'))
	text(numcycles, 14, 'A', col='green')
	text(numcycles, 12, 'C', col='blue')
	text(numcycles, 10, 'G', col='black')
	text(numcycles, 8, 'T', col='yellow')
	text(numcycles, 6, 'N', col='red')

	dev.off()
}

files<-dir(pattern="809_GCCTGAATTTAC_L001_R2_001_AT_QT.paired_matched.fastq.gz")

for(i in 1:length(files)){
  	print(files[i])
  	illumqc(files[i])
}
