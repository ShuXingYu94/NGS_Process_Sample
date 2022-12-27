# Statistics
cd ${workdir}
touch ${stats_dir}/results.txt

# Reads Count (Raw Data) - with seqkit stats
seqkit stats ${basecall_dir}/*.gz -T -j ${threads} > ${stats_dir}/Raw_read_stats.txt
#awk '{print $1,$4}' ${stats_dir}/Raw_read_stats.txt

# Reads Count (After trmming) - with seqkit stats
seqkit stats ${trimmed_dir}/*.gz -T -j ${threads} > ${stats_dir}/Trimmed_read_stats.txt
#${trimmed_dir}+[\S]*R[1|2]up[\S]*[.gz]
# Trimming rate/ Surviving rate


# Mapped Length (After Trimming) - with samtools depth
mkdir ${stats_dir}/Read_depth
for file in $files
do
samtools depth ${aligned_dir}/$file.bam | awk '($3+0) >= 10 {print $0}' > ${stats_dir}/Read_depth/$file.txt
echo -e "$file\tMapped Length\t"`awk "END{print NR}" ${stats_dir}/Read_depth/$file.txt` >> ${stats_dir}/results.txt
done

# Mapped Reads Coverage (After Trimming/ Genome Size)
        #ls ${aligned_dir}/*.bam > ${stats_dir}/file_names.txt
mkdir ${stats_dir}/Coverage
for file in $files
do
samtools coverage ${aligned_dir}/$file.bam -o ${stats_dir}/Coverage/$file.txt
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
# SNPs Count - with bash script using vcf file
# Average snp depth - Figure with vcftools - R
vcftools --vcf ${stacks_dir}/populations.snps.vcf --site-mean-depth  --temp ${log_dir} --out ${stacks_dir}/mean_depth_stat
#proceed in python/R

# snp across chr - python script
if ! [ -n `which python3` ] ; then
  python3 ${stacks_dir}/chr_fig.py
else
  if ! [ -n `which python` ] ; then
    python ${stacks_dir}/chr_fig.py
  else
    echo "Python is not installed. Please install."
  fi
fi
# Mapped Length/SNP


# Consensus Length/SNP


rm -rf ${stats_dir}/bam2txt