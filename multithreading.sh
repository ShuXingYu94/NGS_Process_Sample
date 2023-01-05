#!/bin/sh

# Required Settings
workdir=/Users/.../work_dir
mapping_db=/Users/.../reference_genome.fa
trimmomatic_dir=/Users/.../Trimmomatic-0.39/trimmomatic-0.39.jar
adapter=${workdir}/MIGadapter.fasta
popmap_dir=${workdir}/popmap.txt

# Read file names
files="sample1 sample2"

# Optional settings
illumina_clip=2:30:10
slidingwindow=4:15
threads=6
leading=20
trailing=20
minlen=50
stacks_r=1.0

# Secondary dir
basecall_dir=${workdir}/BaseCall
trimmed_dir=${workdir}/trimmed
aligned_dir=${workdir}/aligned
stacks_dir=${workdir}/stacks
log_dir=${workdir}/log
stats_dir=${workdir}/statistics

# Delete all the files in advance
rm -rf ${workdir}/aligned
rm -rf ${workdir}/log
rm -rf ${workdir}/stacks
rm -rf ${workdir}/trimmed
rm -rf ${workdir}/statistics
rm ${workdir}/*.py
rm ${workdir}/*.r

rm ${workdir}/popmap.txt

# mkdir in workdirection
mkdir ${trimmed_dir}
mkdir ${aligned_dir}
mkdir ${stacks_dir}
mkdir ${log_dir}
mkdir ${stats_dir}

# popmap.txt
if ! test -f "${popmap_dir}"; then
  for file in ${files}
  do
  echo -e "${file}\tpop1" >> ${popmap_dir}
  done
else
  echo "Popmap file already exists, skip automatically creating popmap.txt."
fi

# Fetch file python & R from github
if command -v wget >/dev/null 2>&1; then
  wget -O ${workdir}/chr_fig.py https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/chr_fig.py
  wget -O ${workdir}/MIGadapter.fasta https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/MIGadapter.fasta
else
  if command -v curl >/dev/null 2>&1; then
    curl -o ${workdir}/chr_fig.py https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/chr_fig.py
    curl -o ${workdir}/MIGadapter.fasta https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/MIGadapter.fasta
  else
    echo "Wget and Curl are not installed. Please install."
    exit 0
  fi
fi

# Trimming
for f in ${basecall_dir}/*R1*.gz
do
  {
  out=${f/${basecall_dir}/${trimmed_dir}}
  java -jar ${trimmomatic_dir} PE -phred33 -threads ${threads} ${f} ${f/R1/R2} ${out/R1/R1p} ${out/R1/R1up} ${out/R1/R2p} ${out/R1/R2up} ILLUMINACLIP:${adapter}:${illumina_clip} SLIDINGWINDOW:${slidingwindow} LEADING:${leading} TRAILING:${trailing} MINLEN:${minlen}
  #Optional: -trimlog ${out/.gz/.log}
  } &
done
wait

# Mapping files
for file in ${files}
do
{
  bwa mem -t ${threads} ${mapping_db} ${trimmed_dir}/${file}*R1p*.gz ${trimmed_dir}/${file}*R2p*.gz | samtools view -b | samtools sort --threads {threads} > ${aligned_dir}/${file}.bam
samtools index ${aligned_dir}/${file}.bam
} &
done
wait

# Stacks
gstacks -I ${aligned_dir} -M ${popmap_dir}  -O ${stacks_dir} -t ${threads}
populations -t ${threads} -P ${stacks_dir} -M ${popmap_dir} --structure --vcf  -r ${stacks_r} -O ${stacks_dir}

#!/bin/sh

# Statistics
# Download python and r scripts
if command -v wget >/dev/null 2>&1; then
  wget -O ${workdir}/r_plot.r https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/r_plot.r
else
  if command -v curl >/dev/null 2>&1; then
    curl -o ${workdir}/r_plot.r https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/r_plot.r
  else
    echo "Wget and Curl are not installed. Please install."
    exit 0
  fi
fi

# Create folder and prerequisites
mkdir ${stats_dir}/Coverage
mkdir ${stats_dir}/Read_depth
mkdir ${stats_dir}/Consensus

# Generate depth file
for file in $files
do
{
  if ! test -f "${stats_dir}/Read_depth/$file.txt"; then
  samtools depth ${aligned_dir}/$file.bam | awk '{print $0}' > ${stats_dir}/Read_depth/$file.txt
else
  echo "Depth file exists. Skip creating Read_depth folder"
fi
} &
done
wait

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
echo -e "file\tR1_reads\tR2_reads\tR1_trimmed_reads\tR2_trimmed_reads\tR1_trim_rate\tR2_trim_rate\tTotal_trim_rate\tmean_coverage\tdepth_mean\tmapped_length\treference_length\tconsensus_length\tconsensus_coverage(Genome)\tconsensus_coverage(Mapped)" > ${stats_dir}/results.txt
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
  MAPPED_bases=`awk 'BEGIN{sum=0}{sum += $5}END{print sum}' ${stats_dir}/Coverage/$file.txt`
  DEPTH_mean=`samtools depth ${aligned_dir}/$file.bam | awk '($3+0) > 0 {print $0}' | awk 'BEGIN{sum=0}{sum += $3}END{print sum/NR}'`
  COVERAGE=`awk 'BEGIN{len=0; ref=0}{len += $5; ref += $3}END{print len/ref}' ${stats_dir}/Coverage/$file.txt`

  # Consensus Coverage (Of Genome Size) - with bash script
  # Consensus Coverage (Of Mapped Length After Trimming)) - with bash script
  CONSENSUS_length=`awk 'BEGIN{sum=0}{sum += $2}END{print sum}' ${stats_dir}/Consensus/consensus.txt`
  CONSENSUS_coverage_G=$(printf "%.5f" $(echo "scale=10;${CONSENSUS_length}/${REF_length}" | bc))
  CONSENSUS_coverage_M=$(printf "%.5f" $(echo "scale=10;${CONSENSUS_length}/${MAPPED_bases}" | bc))
  echo -e "$file\t$R1_reads\t$R2_reads\t$R1_trimmed_reads\t$R2_trimmed_reads\t$R1_trim_rate\t$R2_trim_rate\t$Total_trim_rate\t$COVERAGE\t$DEPTH_mean\t$MAPPED_bases\t$REF_length\t$CONSENSUS_length\t$CONSENSUS_coverage_G\t$CONSENSUS_coverage_M" >> ${stats_dir}/results.txt
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
# Whithout filtering
cd ${stacks_dir}
vcftools --vcf populations.snps.vcf --out mean_depth_stat --site-mean-depth
awk '{$4="";print $0}' mean_depth_stat.ldepth.mean |  awk '$0=NR" "$0' > SNP_Mean_Depth.txt
cd ${workdir}
cp ./stacks/SNP_Mean_Depth.txt tmp.txt
R -f r_plot.r
rm tmp.txt
if test -f "Rectangular-Manhattan.MEAN_DEPTH_depth.jpg"; then
  mv Rectangular-Manhattan.MEAN_DEPTH_depth.jpg ./statistics/SNP_Depth_0.jpg
else
  echo "Rscript was not successfully executed."
fi

# With filtering
cd ${stacks_dir}
cat populations.snps.vcf | grep -E '#' > vcfform.txt
sed -e '/#CHROM/r vcfform.txt' snps_depth_10.txt | sed '1d' > snps_depth_10.vcf
rm vcfform.txt
vcftools --vcf snps_depth_10.vcf --site-mean-depth --out mean_depth_stat_10
awk '{$4="";print $0}' mean_depth_stat_10.ldepth.mean |  awk '$0=NR" "$0' > SNP_Mean_Depth_10.txt
cd ${workdir}
cp ./stacks/SNP_Mean_Depth_10.txt tmp.txt
R -f r_plot.r
rm tmp.txt
if test -f "Rectangular-Manhattan.MEAN_DEPTH_depth.jpg"; then
  mv Rectangular-Manhattan.MEAN_DEPTH_depth.jpg ./statistics/SNP_Depth_10.jpg
else
  echo "Rscript was not successfully executed."
fi

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