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
stacks_r=0.8
SNPs_depth=10
read1_symbol=R1
read2_symbol=R2

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
for f in ${basecall_dir}/*${read1_symbol}*.f*
do
out=${f/${basecall_dir}/${trimmed_dir}}
java -jar ${trimmomatic_dir} PE -phred33 -threads ${threads} ${f} ${f/${read1_symbol}/${read2_symbol}} ${out/${read1_symbol}/${read1_symbol}p} ${out/${read1_symbol}/${read1_symbol}up} ${out/${read1_symbol}/${read2_symbol}p} ${out/${read1_symbol}/${read2_symbol}up} ILLUMINACLIP:${adapter}:${illumina_clip} SLIDINGWINDOW:${slidingwindow} LEADING:${leading} TRAILING:${trailing} MINLEN:${minlen}
#Optional: -trimlog ${out/.gz/.log}
done

# Mapping files - Need more fix in regex
for file in ${files}
do
bwa mem -t ${threads} ${mapping_db} ${trimmed_dir}/${file}*${read1_symbol}p*.f* ${trimmed_dir}/${file}*${read2_symbol}p*.f* | samtools view -b | samtools sort --threads {threads} > ${aligned_dir}/${file}.bam
samtools index ${aligned_dir}/${file}.bam
done

# Stacks
gstacks -I ${aligned_dir} -M ${popmap_dir}  -O ${stacks_dir} -t ${threads}
populations -t ${threads} -P ${stacks_dir} -M ${popmap_dir} --structure --vcf  -r ${stacks_r} -O ${stacks_dir}


# Re sample
#ls -1 ${dir} | grep -E "(${file_name}.*(${read1_symbol}|${read2_symbol}).*\.(fastq|fq))|((${read1_symbol}|${read2_symbol}).*${file_name}.*\.(fastq|fq))|((${read1_symbol}|${read2_symbol}).*\.(fastq|fq).*${file_name})"
