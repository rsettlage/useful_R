plotDispEsts<-function(cds)
{
	plot(
		rowMeans(counts(cds,normalized=TRUE)),
		fitInfo(cds)$perGeneDispEsts,
		pch=".",log="xy")

	xg<-10^seq(-.5,5,length.out=300)
	lines(xg,fitInfo(cds)$dispFun(xg),col="red")
}
