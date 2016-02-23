plotDE<-function(res,p){
	plot(
		res$baseMean,
		res$log2FoldChange,
		log="x",pch=20,cex=0.3,
		col=ifelse(res$padj<p,"red","black")
	)
}
