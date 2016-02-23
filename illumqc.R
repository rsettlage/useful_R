illumqc<-function(filename){

require(ShortRead)
fq<-readFastq(filename)
cyclenum<-max(width(fq))
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
base<-strsplit(filename, "\\.")[[1]][1]
bmp<-paste(base, ".bmp", sep="")
bitmap(bmp, res=300)
plot(means, ylab="Mean quality score", xlab="Cycle Number", ylim=c(0, 42), main=lenfq)
arrows(seq(1,max(widfq)), means+sds, seq(1, max(widfq)), means-sds, length=0.01, angle=90, code=3)
abline(h=30, col="red", lty=2)
text(2, 31, labels="High quality",  pos=4)
dev.off()

#Nucleotide profiler

seqmat<-alphabetByCycle(seq)
seqmat2<-rbind(seqmat[1:4,], seqmat[15,])
seqmat3<-seqmat2/lenfq
bmp2<-paste(base, ".np.bmp", collapse="")
bmp2<-sub("\\s+", "", bmp2)
bitmap(bmp2, res=300)
plot(1:cyclenum,seqmat3[1,], type="l", col="green", ylim=c(0, 0.7), ylab="Fraction at each cycle", xlab="Cycle")
lines(1:cyclenum, seqmat3[2,], type="l", col="blue")
lines(1:cyclenum, seqmat3[3,], type="l", col="black")
lines(1:cyclenum, seqmat3[4,], type="l", col="yellow")
lines(1:cyclenum, seqmat3[5,], type="l", col="red")
dev.off()
}

files<-dir(pattern=".gz")

for(i in 1:length(files)){
write(files[i], stdout())
illumqc(files[i])}