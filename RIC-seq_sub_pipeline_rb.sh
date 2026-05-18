mamba activate kingfisher

kingfisher get -r SRR8632820 -m ena-ascp \
 --ascp-ssh-key ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
 --output-directory /home/bioinfo/07_people/rbli1 \
 --download-threads 8 
# error
# force to use key in ascp  # work

conda activate qc_align

BASEpath="/home/bioinfo/07_people/rbli1/RICseq"
SRRid="SRR8632820"
mkdir -p "$BASEpath/Fastqc"
fastqc -o "$BASEpath/Fastqc" \
 --noextract \
 --threads 4 \
 "$BASEpath/FastqFiles/${SRRid}_1.fastq.gz" "$BASEpath/FastqFiles/${SRRid}_2.fastq.gz"  #running!!!!!!!
#添加多线程加速

# 按照去接头-去重-去polyG/N的顺序，因为
# 去poly(N)会改变序列：Cutadapt切掉3'端的N尾巴后，序列变短了，原本是PCR重复的两条reads可能因为切掉的长度不同而不再完全相同
# 创建必要的目录
mkdir -p "$BASEpath/CleanData"
mkdir -p "$BASEpath/Fastqc"

# 步骤1：Trimmomatic 去接头 + 基础质控
# 查找 adapters 文件夹
ls -la $(dirname $(which trimmomatic))/../share/trimmomatic/adapters/
ADAPTER_PATH="/mnt/data/teacher/miniforge3/envs/qc_align/share/trimmomatic/adapters/TruSeq3-PE.fa"

# 步骤1：Trimmomatic - 只去接头，不做质量修剪和长度过滤
trimmomatic PE -threads 4 \
    "$BASEpath/FastqFiles/${SRRid}_1.fastq.gz" \
    "$BASEpath/FastqFiles/${SRRid}_2.fastq.gz" \
    "$BASEpath/CleanData/${SRRid}_trim_R1.fq.gz" \
    "$BASEpath/CleanData/${SRRid}_trim_U1.fq.gz" \
    "$BASEpath/CleanData/${SRRid}_trim_R2.fq.gz" \
    "$BASEpath/CleanData/${SRRid}_trim_U2.fq.gz" \
    ILLUMINACLIP:${ADAPTER_PATH}:2:30:10

# 步骤2：去除 PCR 重复（基于完整序列）
# 把时间信息保存到文件
mkdir -p "$BASEpath/Logs"
( /usr/bin/time -v perl "$BASEpath/remove_duplicated_reads_ricPaper.pl" \
    "$BASEpath/CleanData/${SRRid}_trim_R1.fq.gz" \
    "$BASEpath/CleanData/${SRRid}_trim_R2.fq.gz" \
    "$BASEpath/CleanData/${SRRid}_dedup_R1.fq.gz" \
    "$BASEpath/CleanData/${SRRid}_dedup_R2.fq.gz" ) 2> "$BASEpath/Logs/${SRRid}_dedup_time.log"

# 步骤3：fastp 的 --trim_poly_x 可以处理各种尾巴
fastp -i "$BASEpath/CleanData/${SRRid}_dedup_R1.fq.gz" \    #---------------- here!!!!!!!!
      -I "$BASEpath/CleanData/${SRRid}_dedup_R2.fq.gz" \
      -o "$BASEpath/CleanData/${SRRid}_1_clean.fq.gz" \
      -O "$BASEpath/CleanData/${SRRid}_2_clean.fq.gz" \
      --trim_poly_g \
      --trim_poly_x \
      --length_required 36 \
      --thread 4

# 步骤3：Cutadapt - 处理 polyG + polyN + 长度过滤
cutadapt \
    -a "G{10}" -A "G{10}" \
    -a "N{10}" -A "N{10}" \
    -e 0.1 \
    --minimum-length 36 \
    -o "$BASEpath/CleanData/${SRRid}_1_clean.fq.gz" \
    -p "$BASEpath/CleanData/${SRRid}_2_clean.fq.gz" \
    "$BASEpath/CleanData/${SRRid}_dedup_R1.fq.gz" \
    "$BASEpath/CleanData/${SRRid}_dedup_R2.fq.gz"
#处理 polyG (NovaSeq 平台)
#处理 polyN (RIC-seq 文献要求)
#允许 10% 错配
#**在这里统一做长度过滤**

# 步骤4：FastQC 检查最终结果
fastqc -o "$BASEpath/Fastqc" \
    --noextract \
    --threads 4 \
    "$BASEpath/CleanData/${SRRid}_1_clean.fq.gz" \
    "$BASEpath/CleanData/${SRRid}_2_clean.fq.gz"

# 构建 star-index
conda activate qc_align
conda install -c bioconda star  # 2.7.11b
conda install -c bioconda bowtie2 # version 2.5.5


# 构建 45 pre-rrna bowtie2 索引（.fasta is NR_046235.3）
# bowtie2 快速、内存高效的端到端比对，DNA-seq、ChIP-seq、rRNA 过滤，无内含子，不需要识别剪接位点
#  rRNA 没有内含子
# STAR 处理剪接（跨越内含子）和嵌合体比对，RNA-seq、RIC-seq（需要检测剪接和嵌合）
rRNAgenomPATH="/home/bioinfo/07_people/rbli1/Reference_Genome/rRNA"
mkdir -p "$rRNAgenomPATH/human_45S_pre_rRNA_index_bowtie2"

# bowtie2-build 没有 -o 或 --output 参数，索引文件默认生成在当前目录
cd "$rRNAgenomPATH/human_45S_pre_rRNA_index_bowtie2"
bowtie2-build $rRNAgenomPATH/human_45S_pre_ribosomal_N5.fasta \
 human_45S_pre_rRNA_index
# 过滤 45 pre-rrna
BASEpath="/home/bioinfo/07_people/rbli1/RICseq"
SRRid="SRR8632820"
bowtie2 -x "$rRNAgenomPATH/human_45S_pre_rRNA_index_bowtie2/human_45S_pre_rRNA_index" \
        -1 "$BASEpath/CleanData/${SRRid}_1_clean.fq.gz" \
        -2 "$BASEpath/CleanData/${SRRid}_2_clean.fq.gz" \
        --un-conc-gz "$BASEpath/CleanData/${SRRid}_NOrRNA" \
        --al-conc-gz "$BASEpath/CleanData/${SRRid}_rRNA" \
        --no-unal \
        -p 8 \
        -S /dev/null \
        2> "$BASEpath/${SRRid}_bowtie2.log"
# -x, --un-conc-gz, --al-conc-gz all need provide 前缀

# build star index
genomPATH="/home/bioinfo/07_people/rbli1/Reference_Genome/Human"
mkdir -p "$genomPATH/hg38_index_star2711"
STAR --version # 2.7.11b
STAR --runMode genomeGenerate \
     --genomeDir "$genomPATH/hg38_index_star2711" \
     --genomeFastaFiles "$genomPATH/hg38.fa" \
     --sjdbGTFfile "$genomPATH/hg38.gtf" \
     --runThreadN 8 \
     --genomeSAindexNbases 14   

# Sequence alignment (Chimeric)
# RIC-seq 的 paired-end reads 可能来源于不同 RNA 或跨 RNA–RNA 相互作用的片段，
# 如果直接按普通 paired-end 比对，STAR 会强制配对导致大量 reads 丢失，
# 因此需要 将两端单独比对以保留所有 gapped/chimeric reads。
BASEpath="/home/bioinfo/07_people/rbli1/RICseq"
SRRid="SRR8632820"
genomPATH="/home/bioinfo/07_people/rbli1/Reference_Genome/Human"
mkdir -p "$BASEpath/STARoutput"
STAR --runMode alignReads \
     --genomeDir "$genomPATH/hg38_index_star2711" \
     --readFilesIn "$BASEpath/CleanData/${SRRid}_NOrRNA.1" "$BASEpath/CleanData/${SRRid}_NOrRNA.2" \
     --readFilesCommand zcat \
     --outFileNamePrefix "$BASEpath/STARoutput/${SRRid}_" \
     --outSAMattributes All \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterMultimapNmax 100 \
     --alignIntronMin 1 \
     --scoreGapNoncan -4 \
     --scoreGapATAC -4 \
     --chimSegmentMin 15 \
     --chimJunctionOverhangMin 15 \
     --alignSJoverhangMin 15 \
     --alignSJDBoverhangMin 10 \
     --alignSJstitchMismatchNmax 5 -1 5 5 \
     --chimOutType Junctions WithinBAM SoftClip \
     --runThreadN 16  # running!!!!!

# R1
STAR --runMode alignReads \
     --genomeDir "$genomPATH/hg38_index_star2711" \
     --readFilesIn "$BASEpath/CleanData/${SRRid}_NOrRNA.1" \
     --readFilesCommand zcat \
     --outFileNamePrefix "$BASEpath/STARoutput/${SRRid}_R1_" \
     --outSAMattributes All \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterMultimapNmax 100 \
     --alignIntronMin 1 \
     --scoreGapNoncan -4 \
     --scoreGapATAC -4 \
     --chimSegmentMin 15 \
     --chimJunctionOverhangMin 15 \
     --alignSJoverhangMin 15 \
     --alignSJDBoverhangMin 10 \
     --alignSJstitchMismatchNmax 5 -1 5 5 \
     --chimOutType Junctions SeparateSAMold \
     --runThreadN 16

# R2
STAR --runMode alignReads \
     --genomeDir "$genomPATH/hg38_index_star2711" \
     --readFilesIn "$BASEpath/CleanData/${SRRid}_NOrRNA.2" \
     --readFilesCommand zcat \
     --outFileNamePrefix "$BASEpath/STARoutput/${SRRid}_R2_" \
     --outSAMattributes All \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterMultimapNmax 100 \
     --alignIntronMin 1 \
     --scoreGapNoncan -4 \
     --scoreGapATAC -4 \
     --chimSegmentMin 15 \
     --chimJunctionOverhangMin 15 \
     --alignSJoverhangMin 15 \
     --alignSJDBoverhangMin 10 \
     --alignSJstitchMismatchNmax 5 -1 5 5 \
     --chimOutType Junctions SeparateSAMold \
     --runThreadN 16  #here!!!!!

## Download process code from github of RIC-seq
# https://github.com/caochch/RIC-seq/tree/master/1.find_pairTags_from_STAR_results
# run codes
