# Statistics
cd ${workdir}
touch ${stats_dir}/results.txt

# Reads Count (Raw Data) - with seqkit stats
cd  ${basecall_dir}
echo "#Raw Reads Count" > ${stats_dir}/results.txt
seqkit stats *.gz -T -j ${threads} >> ${stats_dir}/results.txt

# Reads Count (After trmming) - with seqkit stats
cd  ${trimmed_dir}
echo "#Timmed Reads Count" >> ${stats_dir}/results.txt
seqkit stats *R[1-2]p*.gz -T -j ${threads} >> ${stats_dir}/results.txt
cd  ${workdir}

# Trimming rate/ Surviving rate


# Mapped Length (After Trimming) - with samtools depth
#echo "#Mapped Length" >> ${stats_dir}/results.txt
#mkdir ${stats_dir}/Read_depth_0
#for file in $files
#do
#samtools depth ${aligned_dir}/$file.bam | awk '($3+0) > 0 {print $0}' > ${stats_dir}/Read_depth_0/$file.txt
#echo -e "$file\tMapped Length\t"`awk "END{print NR}" ${stats_dir}/Read_depth/$file.txt` >> ${stats_dir}/results.txt
#awk 'BEGIN{sum=0}{sum += $3}END{print sum/NR}' ${stats_dir}/Read_depth/$file.txt >> ${stats_dir}/results.txt
#done

#echo "#Mapped Length with depth >= 10" >> ${stats_dir}/results.txt
#mkdir ${stats_dir}/Read_depth_10
#for file in $files
#do
#samtools depth ${aligned_dir}/$file.bam | awk '($3+0) >= 10 {print $0}' > ${stats_dir}/Read_depth/$file.txt
#echo -e "$file\tMapped Length\t"`awk "END{print NR}" ${stats_dir}/Read_depth/$file.txt` >> ${stats_dir}/results.txt
#done

# Mapped Reads Coverage (After Trimming/ Genome Size)
        #ls ${aligned_dir}/*.bam > ${stats_dir}/file_names.txt
mkdir ${stats_dir}/Coverage
echo "#Coverage" >> ${stats_dir}/results.txt
echo -e "file\tmean_coverage\tdepth_mean\tread_count\tmapped_length\treference_length" >> ${stats_dir}/results.txt
for file in $files
do
samtools coverage ${aligned_dir}/$file.bam -o ${stats_dir}/Coverage/$file.txt
REF_length=`awk 'BEGIN{sum=0}{sum += $3}END{print sum}' ${stats_dir}/Coverage/$file.txt`
READ_sum=`awk 'BEGIN{sum=0}{sum += $4}END{print sum}' ${stats_dir}/Coverage/$file.txt`
MAPPED_bases=`awk 'BEGIN{sum=0}{sum += $5}END{print sum}' ${stats_dir}/Coverage/$file.txt`
DEPTH_mean=`awk 'BEGIN{sum=0}{sum += $3}END{print sum/NR}' ${stats_dir}/Read_depth/$file.txt`
COVERAGE=`awk 'BEGIN{len=0; ref=0}{len += $5; ref += $3}END{print len/ref}' ${stats_dir}/Coverage/$file.txt`
echo -e "$file\t$COVERAGE\t$DEPTH_mean\t$READ_sum\t$MAPPED_bases\t$REF_length" >> ${stats_dir}/results.txt
done

# Average Mapped Reads depth - with samtools depth
# From above file

# Consensus Length - python script
#   awk '{print $1 $2}' ${stats_dir}/Read_depth/d1-LE.txt > ${stats_dir}/Read_depth/d1-LE_ref.txt
#      #awk '{print $1 $2}' ${stats_dir}/Read_depth/d2-LE.txt > ${stats_dir}/Read_depth/d2-LE_ref.txt
#      #grep -wf ${stats_dir}/Read_depth/d1-LE_ref.txt ${stats_dir}/Read_depth/d2-LE_ref.txt
ls ${aligned_dir}/*.bam | awk -F "/" '{print $NF}' | awk -F "." '{print $1}'> ${stats_dir}/file_names.txt

mkdir ${stats_dir}/bam2txt
for file in $files
do
if ! test -f "${aligned_dir}/${file}.txt"; then
  samtools view -@ ${threads} ${aligned_dir}/${file}.bam --threads 6 > ${stats_dir}/bam2txt/${file}.txt
fi
done

if command -v python3 >/dev/null 2>&1; then
  python3 ${workdir}/consensus.py
else
  if command -v python >/dev/null 2>&1; then
    python ${workdir}/consensus.py
  else
    echo "Python is not installed. Please install."
  fi
fi

# Consensus Coverage (Of Genome Size) - with bash script
# Consensus Coverage (Of Mapped Length After Trimming)) - with bash script

echo "#Consensus info" >> ${stats_dir}/results.txt
echo -e "file\tconsensus_length\tconsensus_coverage(Genome)\tconsensus_coverage(Mapped)" >> ${stats_dir}/results.txt
for file in $files
do
CONSENSUS_length=`sed -n '/'$file'/p' ${stats_dir}/compare_output.csv | awk -F "," 'BEGIN{sum=0}{sum += $8}END{print sum}'`
Mapped_length=`sed -n '/'$file'/p' ${stats_dir}/compare_output.csv | awk -F "," 'BEGIN{sum=0}{sum +=$5}END{print sum}'`
REF_length=`awk 'BEGIN{sum=0}{sum += $3}END{print sum}' ${stats_dir}/Coverage/$file.txt`
CONSENSUS_coverage_G=$(printf "%.5f" $(echo "scale=10;${CONSENSUS_length}/${REF_length}" | bc))
CONSENSUS_coverage_M=$(printf "%.5f" $(echo "scale=10;${CONSENSUS_length}/${Mapped_length}" | bc))
echo -e "$file\t$CONSENSUS_length\t$CONSENSUS_coverage_G\t$CONSENSUS_coverage_M" >> ${stats_dir}/results.txt
done

# SNPs Count - with bash script using vcf file

awk -v FS=":" '$8 > 0{ print $0 }' ${workdir}/stacks/populations.snps.vcf | echo -e "SNPs_all\t"`awk "END{print NR}"` >> ${stats_dir}/results.txt
awk -v FS=":" '$8 > 10{ print $0 }' ${workdir}/stacks/populations.snps.vcf | awk -v FS=":" '$12 > 10{ print $0 }' | echo -e "SNPs_depth>10\t"`awk "END{print NR}"` >> ${stats_dir}/results.txt
awk -v FS=":" '$8 > 10{ print $0 }' ${workdir}/stacks/populations.snps.vcf | awk -v FS=":" '$12 > 10{ print $0 }'>  ${stacks_dir}/snps_depth_10.txt

# Mapped Length/SNP
# Consensus Length/SNP
echo "#SNP count" >> ${stats_dir}/results.txt
echo -e "SNP_depth\tSNP_count\tgenome_length_per_SNP\tconsensus_length_per_SNP" >> ${stats_dir}/results.txt

REF_length=`awk 'BEGIN{sum=0}{sum += $3}END{print sum}' ${stats_dir}/Coverage/$file.txt`
CONSENSUS_length=`sed -n '/'$file'/p' ${stats_dir}/compare_output.csv | awk -F "," 'BEGIN{sum=0}{sum += $8}END{print sum}'`

SNP=$(awk -v FS=":" '$8 > 0{ print $0 }' ${workdir}/stacks/populations.snps.vcf | awk "END{print NR}")
SNP_Genome=$(printf "%.1f" $(echo "scale=10;${REF_length}/${SNP}" | bc))
SNP_Consensus=$(printf "%.1f" $(echo "scale=10;${CONSENSUS_length}/${SNP}" | bc))
echo -e "0\t$SNP\t$SNP_Genome\t$SNP_Consensus" >> ${stats_dir}/results.txt

SNP=$(awk -v FS=":" '$8 > 10{ print $0 }' ${workdir}/stacks/populations.snps.vcf | awk "END{print NR}")
SNP_Genome=$(printf "%.1f" $(echo "scale=10;${REF_length}/${SNP}" | bc))
SNP_Consensus=$(printf "%.1f" $(echo "scale=10;${CONSENSUS_length}/${SNP}" | bc))
echo -e "10\t$SNP\t$SNP_Genome\t$SNP_Consensus" >> ${stats_dir}/results.txt

# Average snp depth - Figure with vcftools - R
#vcftools --vcf ${stacks_dir}/populations.snps.vcf --out ${stats_dir}/mean_depth_stat --site-mean-depth
#proceed in python/R

# snp across chr - python script
if command -v python3 >/dev/null 2>&1; then
  python3 ${stacks_dir}/chr_fig.py ${stacks_dir}/populations.snps.vcf ${stats_dir}/SNPs_Distribution_0
  python3 ${stacks_dir}/chr_fig.py ${stacks_dir}/snps_depth_10.txt ${stats_dir}/SNPs_Distribution_10
else
  if command -v python >/dev/null 2>&1; then
    python ${stacks_dir}/chr_fig.py ${stacks_dir}/populations.snps.vcf ${stats_dir}/SNPs_Distribution_0
    python ${stacks_dir}/chr_fig.py ${stacks_dir}/snps_depth_10.txt ${stats_dir}/SNPs_Distribution_10
  else
    echo "Python is not installed. Please install."
  fi
fi



rm -rf ${stats_dir}/bam2txt