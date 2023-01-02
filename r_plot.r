#depth＞0
cd ${stacks_dir}
vcftools --vcf populations.snps.vcf --site-mean-depth | awk '{$4="";print $0}' mean_depth_stat.ldepth.mean |  awk '$0=NR" "$0' >> SNP_Mean_Depth.txt
r -f test.r
if ! test -f "${work_dir}/Rectangular-Manhattan.MEAN_DEPTH_depth.jpg"; then
mv *depth.jpg ${statistics}/depth0.jpg
else `echo rscript was not successfully executed.`
fi

library(CMplot)
data<-read.table(file = "1.txt",header = TRUE)
CMplot(data,plot.type="m",cex=0.8, band=0.5,col=c("grey30","grey60"),threshold = 10,threshold.lty = 1,LOG10=FALSE,ylab = "Average Depth of SNP",file="jpg",memo="depth",dpi=300,file.output=TRUE,verbose=TRUE)

#depth＞10
cat populations.snps.vcf | grep -E '#' > vcfform.txt 
sed -i '1r vcfform.txt' snps_depth_10.txt
sed -i '1d' snps_depth_10.txt
cp snps_depth_10.txt snps_depth_10.vcf
rm vcfform.txt

vcftools --vcf snps_depth_10.vcf --out ./mean_depth_10_stat --site-mean-depth
awk '{$4="";print $0}' mean_depth_10_stat.ldepth.mean |  awk '$0=NR" "$0' >> 2.txt



library(CMplot)
data<-read.table(file = "2.txt",header = TRUE)
CMplot(data,plot.type="m",cex=0.8, band=0.5,col=c("grey30","grey60"),threshold = 10,threshold.lty = 1,LOG10=FALSE,ylab = "Average Depth of SNP",file="jpg",memo="depth10",dpi=300,file.output=TRUE,verbose=TRUE)
