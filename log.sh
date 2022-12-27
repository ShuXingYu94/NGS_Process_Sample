#!/bin/bash

# Required Settings
workdir=/Users/zhuxingyu/20221218_MIG_dpMIG/test_dir
mapping_db=/Users/zhuxingyu/Reference_Genome/BnapusDarmor-bzh/BnapusDarmor-bzh.fa
trimmomatic_dir=/Users/zhuxingyu/ADS/Trimmomatic-0.39/trimmomatic-0.39.jar
adapter=${workdir}/MIGadapter.fasta
popmap_dir=${workdir}/popmap.txt

# Read file names
files="d1-LE
d2-LE"

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

# mkdir in workdirection
mkdir ${trimmed_dir}
mkdir ${aligned_dir}
mkdir ${stacks_dir}
mkdir ${log_dir}
mkdir ${stats_dir}

# popmap.txt
if [ ! -d ${popmap_dir} ]; then
  for file in ${files}
  do
  echo -e "${file}\tpop1" >> ${popmap_dir}
  done
fi

# Fetch file python & R from github
if ! [ -n `which wget` ] ; then
  wget -O ${workdir}/chr_fig.py https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/chr_fig.py
  wget -O ${workdir}/consensus.py https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/consensus.py
else
  if [ -n `which curl` ] ; then
    curl -o ${workdir}/chr_fig.py https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/chr_fig.py
    curl -o ${workdir}/consensus.py https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/consensus.py
  else
    echo "Wget and Curl are not installed. Please install."
    exit 0
  fi
fi

# Trimming
for f in ${basecall_dir}/*R1*.gz
do
out=${f/${basecall_dir}/${trimmed_dir}}
java -jar ${trimmomatic_dir} PE -phred33 -threads ${threads} ${f} ${f/R1/R2} ${out/R1/R1p} ${out/R1/R1up} ${out/R1/R2p} ${out/R1/R2up} ILLUMINACLIP:${adapter}:${illumina_clip} SLIDINGWINDOW:${slidingwindow} LEADING:${leading} TRAILING:${trailing} MINLEN:${minlen}
#Optional: -trimlog ${out/.gz/.log}
done

# Mapping files
for file in ${files}
do
bwa mem -t ${threads} ${mapping_db} ${trimmed_dir}/${file}*R1p*.gz ${trimmed_dir}/${file}*R2p*.gz | samtools view -b | samtools sort --threads {threads} > ${aligned_dir}/${file}.bam
samtools index ${aligned_dir}/${file}.bam
done

# Stacks
gstacks -I ${aligned_dir} -M ${popmap_dir}  -O ${stacks_dir} -t ${threads}
populations -t ${threads} -P ${stacks_dir} -M ${popmap_dir} --structure --vcf  -r ${stacks_r} -O ${stacks_dir}

# Statistics

# Reads Count (Raw Data) - with seqkit stats
seqkit stats ${basecall_dir}/*.gz -T -j ${threads} > ${stats_dir}/Raw_read_stats.txt
#awk '{print $1,$4}' ${stats_dir}/Raw_read_stats.txt

# Reads Count (After trmming) - with seqkit stats
seqkit stats ${trimmed_dir}/*.gz -T -j ${threads} > ${stats_dir}/Trimmed_read_stats.txt
#${trimmed_dir}+[\S]*R[1|2]up[\S]*[.gz]
# Trimming rate/ Surviving rate

# Mapped Length (After Trimming) - with samtools
samtools view test.bam | awk '{print length($10)}'
samtools view ${aligned_dir}/d1-LE.bam | awk '{print length($10)}'
# Mapped Reads Coverage (After Trimming/ Genome Size)
#   Need Genome Size(from chr file?) - with samtools
# Average Mapped Reads depth - with samtools

# Consensus Length - python script
ls ${aligned_dir}/*.txt> ${stats_dir}/file_names.txt

for file in $files
do
if ! test -f "${stats_dir}/${file}.txt"; then
  samtools view -@ ${threads} ${aligned_dir}/${file}.bam --threads 6 > ${stats_dir}/${file}.txt
fi
done

# Consensus Coverage (Of Genome Size) - with bash script
# Consensus Coverage (Of Mapped Length After Trimming)) - with bash script
# SNPs Count - with bash script using vcf file
# Average snp depth - Figure with vcftools - R
vcftools --vcf ${stacks_dir}/populations.snps.vcf --site-mean-depth  --temp ${log_dir} --out ${stacks_dir}/mean_depth_stat
#proceed in python

# snp across chr - python script
python3 ${stacks_dir}/chr_fig.py
# Mapped Length/SNP


# Consensus Length/SNP
