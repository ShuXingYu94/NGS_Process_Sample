#!/bin/sh

# Settings
draw_genome_border="FALSE" # FALSE if genoome border is not needed in distribution figure.

# Statistics
export SNPs_depth

# Download python and r scripts
cd ${workdir}
if command -v wget >/dev/null 2>&1; then
  wget -O ${workdir}/r_plot.r https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/r_plot.r
  wget -O ${workdir}/r_distribution.r https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/r_distribution.r

else
  if command -v curl >/dev/null 2>&1; then
    curl -o ${workdir}/r_plot.r https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/r_plot.r
    curl -o ${workdir}/r_distribution.r https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/r_distribution.r
  else
    echo "Wget and Curl are not installed. Please install."
    exit 0
  fi
fi

# Create folder and prerequisites
if ! test -f "${stats_dir}"; then
  mkdir ${stats_dir}
fi
touch ${stats_dir}/results.txt
mkdir ${stats_dir}/Coverage
mkdir ${stats_dir}/Read_depth
mkdir ${stats_dir}/Consensus

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
num_files=`echo $files | wc -w`
filter_depth=`echo "scale=0; $num_files * $stacks_r" | bc`
export filter_depth

if ! test -f "${stats_dir}/Consensus/consensus.txt"; then
  for chr in $chromosomes
  do
    {
    cat ${stats_dir}/Read_depth/*.txt | grep $chr | awk '{print $2}' | sort | uniq -c > ${log_dir}/consensus_${chr}.txt
    CONSENSUS_length=`awk '$1 >= ENVIRON["filter_depth"]{ print $0 }' ${log_dir}/consensus_${chr}.txt | awk "END{print NR}"`
    map_length=`awk 'END{print NR}' ${log_dir}/consensus_${chr}.txt`
    echo "$chr" "${CONSENSUS_length}" "$map_length" >> ${stats_dir}/Consensus/tmp_consensus.txt
    } &
  done
  wait
  sort ${stats_dir}/Consensus/tmp_consensus.txt > ${stats_dir}/Consensus/consensus.txt
  rm ${stats_dir}/Consensus/tmp_consensus.txt
  rm ${log_dir}/consensus_*.txt
else
  echo "Consensus sequence file exists. Skip creating consensus length calculation."
fi

# Main process
echo -e "file\tR1_reads\tR2_reads\tR1_trimmed_reads\tR2_trimmed_reads\tR1_trim_rate\tR2_trim_rate\tTotal_trim_rate\tmean_coverage\tdepth_mean\tmapped_length\total_mapped_length\treference_length\tconsensus_length\tconsensus_coverage(Genome)\tconsensus_coverage(Mapped)" > ${stats_dir}/results.txt
for file in $files
do
  echo "Analyzing file ${file}."

  # Reads Count (Raw Data) - with seqkit stats
  cd ${basecall_dir}
  fn1=`ls -1 ${basecall_dir} | grep -E "(.*${file}.*${read1_symbol}.*\.(fastq|fq))|(.*${read1_symbol}.*${file}.*\.(fastq|fq))|(.*${read1_symbol}.*\.(fastq|fq).*${file})"`
  fn2=`ls -1 ${basecall_dir} | grep -E "(.*${file}.*${read2_symbol}.*\.(fastq|fq))|(.*${read2_symbol}.*${file}.*\.(fastq|fq))|(.*${read2_symbol}.*\.(fastq|fq).*${file})"`
  if ! `$fn1 | awk '{print NR}'`==1; then
    echo "ERROR: multiple (or none) file name containing ${file}, ${read1_symbol} and fastq or fq are(is) found in directory '${basecall_dir}'."
  fi
  R1_reads=`seqkit stats $fn1 -T -j ${threads} | awk 'NR>1{print $4}'`
  R2_reads=`seqkit stats $fn2 -T -j ${threads} | awk 'NR>1{print $4}'`

  echo "Reads Count (Raw Data) parsed."

  # Reads Count (After trmming) - with seqkit stats
  cd  ${trimmed_dir}
  fn3=`ls -1 ${trimmed_dir} | grep -E "(.*${file}.*${read1_symbol}p.*\.(fastq|fq))|(.*${read1_symbol}p.*${file}.*\.(fastq|fq))|(.*${read1_symbol}p.*\.(fastq|fq).*${file})"`
  fn4=`ls -1 ${trimmed_dir} | grep -E "(.*${file}.*${read2_symbol}p.*\.(fastq|fq))|(.*${read2_symbol}p.*${file}.*\.(fastq|fq))|(.*${read2_symbol}p.*\.(fastq|fq).*${file})"`
  if ! `$fn3 | awk '{print NR}'`==1; then
    echo "ERROR: multiple(or none) file name containing '${file}', '${read1_symbol}p' and 'fastq' or 'fq' are(is) found in directory '${trimmed_dir}'."
  fi
  R1_trimmed_reads=`seqkit stats $fn3 -T -j ${threads} | awk 'NR>1{print $4}'`
  R2_trimmed_reads=`seqkit stats $fn4 -T -j ${threads} | awk 'NR>1{print $4}'`
  R1_trim_rate=$(printf "%.5f" $(echo "scale=10;${R1_trimmed_reads}/${R1_reads}" | bc))
  R2_trim_rate=$(printf "%.5f" $(echo "scale=10;${R2_trimmed_reads}/${R2_reads}" | bc))
  Total_trim_rate=$(printf "%.5f" $(echo "scale=10;(${R1_trimmed_reads}+${R2_trimmed_reads})/(${R1_reads}+${R2_reads})" | bc))
  cd  ${workdir}

  echo "Reads Count (After trmming) parsed."

  # Mapped Reads Coverage (After Trimming/ Genome Size)
  if ! test -f "${stats_dir}/Coverage/$file.txt"; then
    samtools coverage ${aligned_dir}/$file.bam -o ${stats_dir}/Coverage/$file.txt
  else
    echo "Coverage file exists. Skip creating coverage calculation."
  fi
  REF_length=`awk 'BEGIN{sum=0}{sum += $3}END{print sum}' ${stats_dir}/Coverage/$file.txt`
  MAPPED_bases=`awk 'BEGIN{sum=0}{sum += $5}END{print sum}' ${stats_dir}/Coverage/$file.txt`
  DEPTH_mean=`samtools depth ${aligned_dir}/$file.bam | awk '($3+0) > 0 {print $0}' | awk 'BEGIN{sum=0}{sum += $3}END{print sum/NR}'`
  COVERAGE=`awk 'BEGIN{len=0; ref=0}{len += $5; ref += $3}END{print len/ref}' ${stats_dir}/Coverage/$file.txt`

  echo "Mapped Reads Coverage data parsed."

  # Consensus Coverage (Of Genome Size) - with bash script
  # Consensus Coverage (Of Mapped Length After Trimming)) - with bash script
  CONSENSUS_length=`awk 'BEGIN{sum=0}{sum += $2}END{print sum}' ${stats_dir}/Consensus/consensus.txt`
  total_mapped_length=`awk 'BEGIN{sum=0}{sum += $3}END{print sum}' ${stats_dir}/Consensus/consensus.txt`
  CONSENSUS_coverage_G=$(printf "%.5f" $(echo "scale=10;${CONSENSUS_length}/${REF_length}" | bc))
  CONSENSUS_coverage_M=$(printf "%.5f" $(echo "scale=10;${CONSENSUS_length}/${total_mapped_length}" | bc))

  echo "Consensus analysis data parsed."

  echo -e "$file\t$R1_reads\t$R2_reads\t$R1_trimmed_reads\t$R2_trimmed_reads\t$R1_trim_rate\t$R2_trim_rate\t$Total_trim_rate\t$COVERAGE\t$DEPTH_mean\t$MAPPED_bases\t$total_mapped_length\t$REF_length\t$CONSENSUS_length\t$CONSENSUS_coverage_G\t$CONSENSUS_coverage_M" >> ${stats_dir}/results.txt
done

# SNPs Count - with bash script using vcf file
# Mapped Length/SNP
# Consensus Length/SNP

awk '/#CHROM/{print $0}' ${workdir}/stacks/populations.snps.vcf > ${stacks_dir}/snps_depth_${SNPs_depth}.txt
awk -v FS=":" '$8 >= ENVIRON["SNPs_depth"]{ print $0 }' ${workdir}/stacks/populations.snps.vcf | awk -v FS=":" '$12 >= ENVIRON["SNPs_depth"]{ print $0 }' >> ${stacks_dir}/snps_depth_${SNPs_depth}.txt

echo -e "\n#SNP count" >> ${stats_dir}/results.txt
echo -e "SNP_depth\tSNP_count\tgenome_length_per_SNP\tconsensus_length_per_SNP" >> ${stats_dir}/results.txt

REF_length=`awk 'BEGIN{sum=0}{sum += $3}END{print sum}' ${stats_dir}/Coverage/$file.txt`
CONSENSUS_length=`awk 'BEGIN{sum=0}{sum += $2}END{print sum}' ${stats_dir}/Consensus/consensus.txt`

SNP=$(awk -v FS=":" '$8 > 0{ print $0 }' ${workdir}/stacks/populations.snps.vcf | awk "END{print NR}")
SNP_Genome=$(printf "%.1f" $(echo "scale=10;${REF_length}/${SNP}" | bc))
SNP_Consensus=$(printf "%.1f" $(echo "scale=10;${CONSENSUS_length}/${SNP}" | bc))
echo -e "SNPs depth: 0\t$SNP\t$SNP_Genome\t$SNP_Consensus" >> ${stats_dir}/results.txt

SNP=$(awk -v FS=":" '$8 >=ENVIRON["SNPs_depth"]{ print $0 }' ${workdir}/stacks/populations.snps.vcf | awk -v FS=":" '$12 >= ENVIRON["SNPs_depth"]{ print $0 }' | awk "END{print NR}")
SNP_Genome=$(printf "%.1f" $(echo "scale=10;${REF_length}/${SNP}" | bc))
SNP_Consensus=$(printf "%.1f" $(echo "scale=10;${CONSENSUS_length}/${SNP}" | bc))
echo -e "SNPs depth: $SNPs_depth\t$SNP\t$SNP_Genome\t$SNP_Consensus" >> ${stats_dir}/results.txt

# Get figures
# Average snp depth - Figure with vcftools - R
# snp across chr - R
# Whithout filtering
cd ${stacks_dir}
vcftools --vcf populations.snps.vcf --out mean_depth_stat --site-mean-depth
awk '{$4="";print $0}' mean_depth_stat.ldepth.mean |  awk '$0=NR" "$0' > SNP_Mean_Depth.txt
cd ${workdir}
cp ./stacks/SNP_Mean_Depth.txt tmp.txt
if ! draw_genome_border="FALSE";then
  awk 'NR>1{print $2,$1,$3,0}' ./statistics/Coverage/d1-LE.txt >> tmp.txt # add the right border
fi
R -f r_plot.r
R -f r_distribution.r
rm tmp.txt
if test -f "Rectangular-Manhattan.MEAN_DEPTH_depth.jpg"; then
  mv Rectangular-Manhattan.MEAN_DEPTH_depth.jpg ./statistics/SNP_Depth_0.jpg
else
  echo "SNP_Depth_0.jpg was not successfully created."
fi
if test -f "SNP-Density.MEAN_DEPTH.jpg"; then
  mv SNP-Density.MEAN_DEPTH.jpg ./statistics/SNPs_Distribution_0.jpg
else
  echo "SNPs_Distribution_0.jpg was not successfully created."
fi


# With filtering
cd ${stacks_dir}
cat populations.snps.vcf | grep -E '#' > vcfform.txt
sed -e '/#CHROM/r vcfform.txt' snps_depth_${SNPs_depth}.txt | sed '1d' > snps_depth_${SNPs_depth}.vcf
rm vcfform.txt
vcftools --vcf snps_depth_${SNPs_depth}.vcf --site-mean-depth --out mean_depth_stat_${SNPs_depth}
awk '{$4="";print $0}' mean_depth_stat_${SNPs_depth}.ldepth.mean |  awk '$0=NR" "$0' > SNP_Mean_Depth_${SNPs_depth}.txt
cd ${workdir}
cp ./stacks/SNP_Mean_Depth_${SNPs_depth}.txt tmp.txt
if ! draw_genome_border="FALSE";then
  awk 'NR>1{print $2,$1,$3,0}' ./statistics/Coverage/d1-LE.txt >> tmp.txt # add the right border
fi
R -f r_plot.r
R -f r_distribution.r
rm tmp.txt
if test -f "Rectangular-Manhattan.MEAN_DEPTH_depth.jpg"; then
  mv Rectangular-Manhattan.MEAN_DEPTH_depth.jpg ./statistics/SNP_Depth_${SNPs_depth}.jpg
else
  echo "SNP_Depth_${SNPs_depth}.jpg was not successfully created."
fi
if test -f "SNP-Density.MEAN_DEPTH.jpg"; then
  mv SNP-Density.MEAN_DEPTH.jpg ./statistics/SNPs_Distribution_${SNPs_depth}.jpg
else
  echo "SNPs_Distribution_${SNPs_depth}.jpg was not successfully created."
fi