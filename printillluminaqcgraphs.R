##################################################################
#
#	Author: K. Wyatt McMahon, Ph.D.
#
#	Date: 12/19/2011
#
#	printilluminaqcgraphs.R - converts.ilqc files to graphs
#
#	usage: R BATCH printqcfiles.R
###################################################################
plot_qvals<-function(vector, filename){
vals<-vector[which(vector>0)]
bpname<-paste(c("qplot", filename, ".bmp"), collapse="_")
bitmap(bpname, res=300)
plot(1:length(vals), vals, ylab="Mean quality score", xlab="Cycle number", ylim=c(0, 40))
lines(1:length(vals), rep(30, length(vals)), lty=2, col="red")
text(0, 29.5, labels="High quality", pos=4, lwd=5) 
dev.off()
}

plot_profile<-function(data, filename){
data<-data[,which(!is.na(data[3,]))]
len<-ncol(data)
totals<-apply(data[3:7,], 2, sum, na.rm=T)
as<-data[3,]/totals
gs<-data[4,]/totals
cs<-data[5,]/totals
ts<-data[6,]/totals
ns<-data[7,]/totals

if(sum(ns, na.rm=T)==0){
ns<-rep(0, len)
}

ns2<-make_profile_vector(data[7,])
if(sum(ns2, na.rm=T)==0){
ns2<-rep(0, len)
}

prname<-paste(c("profile", filename, ".bmp"), collapse="_")
bitmap(prname, res=300)
plot(1:len, as, ylab="Fraction at position", xlab="Position", col=c("green"), type="l", ylim=c(0, 0.5))
lines(1:len, ns, col="red", type="l")
lines(1:len, gs, col="black", type="l")
lines(1:len, cs, col="blue", type="l")
lines(1:len, ts, col="orange", type="l")

lines(1:len, ns2, col="red", type="l", lty=2)
legend(len-(0.25*len), 0.5, legend=c("N", "A", "G", "C", "T", "Ns per total N"), col=c("red", "green", "black", "blue", "orange"), lty=c(rep(1, 5), 2))
dev.off()
}



make_profile_vector<-function(vector){
nas<-is.na(vector)
vector[nas]<-0

sumv<-sum(vector, na.rm=T)
vector<-vector/sumv
return(vector)
}


files<-dir()

files<-files[which(regexpr(".out$", files, perl=T)>0)]


for(i in 1:length(files)){

	data<-read.delim(files[i], header=F, stringsAsFactors=F)
	plot_qvals(data[2,], files[i])
	plot_profile(data, files[i])

}
