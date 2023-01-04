library(CMplot);
data<-read.table(file = "./tmp.txt",header = TRUE);
CMplot(data,plot.type="m",cex=0.8, band=0.5,col=c("grey30","grey60"),threshold = 10,threshold.lty = 1,LOG10=FALSE,ylab = "Average Depth of SNP", file="jpg", memo="depth", dpi=300, file.output=TRUE, verbose=TRUE);