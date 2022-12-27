# Statistics

# Reads Count (Raw Data) - with seqkit stats
seqkit stats ${basecall_dir}/*.gz -T -j ${threads} > ${stats_dir}/Raw_read_stats.txt
#awk '{print $1,$4}' ${stats_dir}/Raw_read_stats.txt

# Reads Count (After trmming) - with seqkit stats
seqkit stats ${trimmed_dir}/*.gz -T -j ${threads} > ${stats_dir}/Trimmed_read_stats.txt
#${trimmed_dir}+[\S]*R[1|2]up[\S]*[.gz]
# Trimming rate/ Surviving rate

# Mapped Length (After Trimming) - with samtools depth
samtools view test.bam | awk '{print length($10)}'
samtools view ${aligned_dir}/d1-LE.bam | awk '{print length($10)}'

# Mapped Reads Coverage (After Trimming/ Genome Size)
ls ${aligned_dir}/*.bam > ${stats_dir}/file_names.txt
samtools coverage -b ${stats_dir}/file_names.txt

#   Need Genome Size(from chr file?) - with samtools
# Average Mapped Reads depth - with samtools depth

# Consensus Length - python script
ls ${aligned_dir}/*.bam > ${stats_dir}/file_names.txt
samtools coverage -b ${stats_dir}/file_names.txt


for file in $files
do
if ! test -f "${aligned_dir}/${file}.txt"; then
  samtools view -@ ${threads} ${aligned_dir}/${file}.bam --threads 6 > ${aligned_dir}/${file}.txt
fi
done

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
  if [ -n `which python` ] ; then
    python ${stacks_dir}/chr_fig.py
  else
    echo "Python is not installed. Please install."
  fi
fi
# Mapped Length/SNP


# Consensus Length/SNP
