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

# mkdir in workdirection
mkdir ${trimmed_dir}
mkdir ${aligned_dir}
mkdir ${stacks_dir}
mkdir ${log_dir}

# popmap.txt
for file in ${files}
do
echo -e "${file}\tpop1">> ${popmap_dir}
done

for file in ${files}
do
echo -e "${file}\tpop1"
done

# Fetch file python & R from github
# Need special from macos or windows/linux systems

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
# Reads Count (After trmming) - with seqkit stats
# Trimming rate/ Surviving rate
# Mapped Length (After Trimming) - with samtools
# Mapped Reads Coverage (After Trimming/ Genome Size)
#   Need Genome Size(from chr file?) - with samtools
# Average Mapped Reads depth - with samtools
# Consensus Length - python script
# Consensus Coverage (Of Genome Size) - with bash script
# Consensus Coverage (Of Mapped Length After Trimming)) - with bash script
# SNPs Count - with bash script using vcf file
# Average snp depth - Figure with vcftools - R
vcftools --vcf ${stacks_dir}/populations.snps.vcf --site-mean-depth  --temp ${log_dir} --out ${stacks_dir}/mean_depth_stat
proceed in python
# snp across chr - python script
python3 ${stacks_dir}/chr_fig.py
# Mapped Length/SNP
# Consensus Length/SNP



# for test
rm -rf ${workdir}/aligned
rm -rf ${workdir}/log
rm ${workdir}/popmap.txt
rm -rf ${workdir}/stacks
rm -rf ${workdir}/trimmed