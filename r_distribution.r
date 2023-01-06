library(CMplot);
data<-read.table(file = "./tmp.txt",header = TRUE);
CMplot(data,type="p",plot.type="d",bin.size=1e3,chr.den.col=c("grey","black","red"),file="jpg",memo="",dpi=300,main="",file.output=TRUE,verbose=TRUE,width=9,height=6)
