#!/bin/sh

# Statistics
# Create folder and prerequisites
touch ${stats_dir}/results.txt
mkdir ${stats_dir}/Coverage
mkdir ${stats_dir}/Read_depth
mkdir ${stats_dir}/Consensus
#touch ${stats_dir}/Consensus/consensus.txt
#touch ${stats_dir}/Consensus/${file}_mapped_length.txt

# Generate depth file
for file in $files
do
if ! test -f "${stats_dir}/Read_depth/$file.txt"; then
  samtools depth ${aligned_dir}/$file.bam | awk '{print $0}' > ${stats_dir}/Read_depth/$file.txt
else
  echo "Depth file exists. Skip creating Read_depth folder"
fi
done

# Get chromosome file
awk -v FS="\t"  'NR> 1{print $1}' ${stacks_dir}/catalog.chrs.tsv > ${stats_dir}/Consensus/chr.txt
chromosomes=`cat ${stats_dir}/Consensus/chr.txt`

# consensus length
if ! test -f "${stats_dir}/Consensus/consensus.txt"; then
  for chr in $chromosomes
  do
    {
      CONSENSUS_length=`cat ${stats_dir}/Read_depth/*.txt | grep $chr | awk '{print $2}' | sort | uniq -cd | awk '$1 > 1{ print $0 }' | awk "END{print NR}"`
      echo "$chr" "${CONSENSUS_length}" >> ${stats_dir}/Consensus/tmp_consensus.txt
      } &
  done
  wait
  sort ${stats_dir}/Consensus/tmp_consensus.txt > ${stats_dir}/Consensus/consensus.txt
  rm ${stats_dir}/Consensus/tmp_consensus.txt
else
  echo "Consensus sequence file exists. Skip creating consensus length calculation."
fi

# Main process
echo -e "file\tR1_reads\tR2_reads\tR1_trimmed_reads\tR2_trimmed_reads\tR1_trim_rate\tR2_trim_rate\tTotal_trim_rate\tmean_coverage\tdepth_mean\tread_count\tmapped_length\treference_length\tconsensus_length\tconsensus_coverage(Genome)\tconsensus_coverage(Mapped)" > ${stats_dir}/results.txt
for file in $files
do
  # Reads Count (Raw Data) - with seqkit stats
  cd ${basecall_dir}
  fn1=`ls ${basecall_dir} | grep *$file*R1*.gz`
  fn2=`ls ${basecall_dir} | grep *$file*R1*.gz`
  R1_reads=`seqkit stats $fn1 -T -j ${threads} | awk 'NR>1{print $4}'`
  R2_reads=`seqkit stats $fn2 -T -j ${threads} | awk 'NR>1{print $4}'`

  # Reads Count (After trmming) - with seqkit stats
  cd  ${trimmed_dir}
  fn3=`ls ${trimmed_dir} | grep $file |grep R1p`
  fn4=`ls ${trimmed_dir} | grep $file |grep R2p`
  #seqkit stats *R[1-2]p*.gz -T -j ${threads} >> ${stats_dir}/tmp2.txt
  R1_trimmed_reads=`seqkit stats $fn3 -T -j ${threads} | awk 'NR>1{print $4}'`
  R2_trimmed_reads=`seqkit stats $fn4 -T -j ${threads} | awk 'NR>1{print $4}'`
  R1_trim_rate=$(printf "%.5f" $(echo "scale=10;${R1_trimmed_reads}/${R1_reads}" | bc))
  R2_trim_rate=$(printf "%.5f" $(echo "scale=10;${R2_trimmed_reads}/${R2_reads}" | bc))
  Total_trim_rate=$(printf "%.5f" $(echo "scale=10;(${R1_trimmed_reads}+${R2_trimmed_reads})/(${R1_reads}+${R2_reads})" | bc))
  cd  ${workdir}

  # Mapped Reads Coverage (After Trimming/ Genome Size)
  if ! test -f "${stats_dir}/Coverage/$file.txt"; then
    samtools coverage ${aligned_dir}/$file.bam -o ${stats_dir}/Coverage/$file.txt
  else
    echo "Coverage file exists. Skip creating coverage calculation."
  fi
  REF_length=`awk 'BEGIN{sum=0}{sum += $3}END{print sum}' ${stats_dir}/Coverage/$file.txt`
  READ_sum=`awk 'BEGIN{sum=0}{sum += $4}END{print sum}' ${stats_dir}/Coverage/$file.txt`
  MAPPED_bases=`awk 'BEGIN{sum=0}{sum += $5}END{print sum}' ${stats_dir}/Coverage/$file.txt`
  DEPTH_mean=`samtools depth ${aligned_dir}/$file.bam | awk '($3+0) > 0 {print $0}' | awk 'BEGIN{sum=0}{sum += $3}END{print sum/NR}'`
  COVERAGE=`awk 'BEGIN{len=0; ref=0}{len += $5; ref += $3}END{print len/ref}' ${stats_dir}/Coverage/$file.txt`

  # Consensus Coverage (Of Genome Size) - with bash script
  # Consensus Coverage (Of Mapped Length After Trimming)) - with bash script
  CONSENSUS_length=`awk 'BEGIN{sum=0}{sum += $2}END{print sum}' ${stats_dir}/Consensus/consensus.txt`
  CONSENSUS_coverage_G=$(printf "%.5f" $(echo "scale=10;${CONSENSUS_length}/${REF_length}" | bc))
  CONSENSUS_coverage_M=$(printf "%.5f" $(echo "scale=10;${CONSENSUS_length}/${MAPPED_bases}" | bc))
  echo -e "$file\t$R1_reads\t$R2_reads\t$R1_trimmed_reads\t$R2_trimmed_reads\t$R1_trim_rate\t$R2_trim_rate\t$Total_trim_rate\t$COVERAGE\t$DEPTH_mean\t$READ_sum\t$MAPPED_bases\t$REF_length\t$CONSENSUS_length\t$CONSENSUS_coverage_G\t$CONSENSUS_coverage_M" >> ${stats_dir}/results.txt
done

# SNPs Count - with bash script using vcf file
# Mapped Length/SNP
# Consensus Length/SNP
awk '/#CHROM/{print $0}' ${workdir}/stacks/populations.snps.vcf > ${stacks_dir}/snps_depth_10.txt
awk -v FS=":" '$8 > 10{ print $0 }' ${workdir}/stacks/populations.snps.vcf | awk -v FS=":" '$12 > 10{ print $0 }' >> ${stacks_dir}/snps_depth_10.txt


echo -e "\n#SNP count" >> ${stats_dir}/results.txt
echo -e "SNP_depth\tSNP_count\tgenome_length_per_SNP\tconsensus_length_per_SNP" >> ${stats_dir}/results.txt

REF_length=`awk 'BEGIN{sum=0}{sum += $3}END{print sum}' ${stats_dir}/Coverage/$file.txt`
CONSENSUS_length=`awk 'BEGIN{sum=0}{sum += $2}END{print sum}' ${stats_dir}/Consensus/consensus.txt`

SNP=$(awk -v FS=":" '$8 > 0{ print $0 }' ${workdir}/stacks/populations.snps.vcf | awk "END{print NR}")
SNP_Genome=$(printf "%.1f" $(echo "scale=10;${REF_length}/${SNP}" | bc))
SNP_Consensus=$(printf "%.1f" $(echo "scale=10;${CONSENSUS_length}/${SNP}" | bc))
echo -e "0\t$SNP\t$SNP_Genome\t$SNP_Consensus" >> ${stats_dir}/results.txt

SNP=$(awk -v FS=":" '$8 > 10{ print $0 }' ${workdir}/stacks/populations.snps.vcf | awk -v FS=":" '$12 > 10{ print $0 }' | awk "END{print NR}")
SNP_Genome=$(printf "%.1f" $(echo "scale=10;${REF_length}/${SNP}" | bc))
SNP_Consensus=$(printf "%.1f" $(echo "scale=10;${CONSENSUS_length}/${SNP}" | bc))
echo -e "10\t$SNP\t$SNP_Genome\t$SNP_Consensus" >> ${stats_dir}/results.txt


# Get figures

# Average snp depth - Figure with vcftools - R
#vcftools --vcf ${stacks_dir}/populations.snps.vcf --out ${stats_dir}/mean_depth_stat --site-mean-depth
#proceed in python/R

# snp across chr - python script
cd ${workdir}
if command -v python3 >/dev/null 2>&1; then
  python3 ${workdir}/chr_fig.py ${stacks_dir}/populations.snps.vcf ${stats_dir}/SNPs_Distribution_0
  python3 ${workdir}/chr_fig.py ${stacks_dir}/snps_depth_10.txt ${stats_dir}/SNPs_Distribution_10
else
  if command -v python >/dev/null 2>&1; then
    python ${workdir}/chr_fig.py ${stacks_dir}/populations.snps.vcf ${stats_dir}/SNPs_Distribution_0
    python ${workdir}/chr_fig.py ${stacks_dir}/snps_depth_10.txt ${stats_dir}/SNPs_Distribution_10
  else
    echo "Python is not installed. Please install."
  fi
fi