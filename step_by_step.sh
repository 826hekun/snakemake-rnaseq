# 单个样本 fastp
mamba install -c bioconda fastp
fastp -i raw_data/Cs1_rep1_1.fastq.gz \
        -I raw_data/Cs1_rep1_2.fastq.gz \
        -o clean_data/Cs1_rep1_1.fastp.fastq.gz \
        -O clean_data/Cs1_rep1_2.fastp.fastq.gz \
        -h clean_data/Cs1_rep1.fastp.html \
        -j clean_data/Cs1_rep1.fastp.json

# 构建参考基因组
mamba install hisat2
hisat2-build ref/genome.fa genome.hisat2

# 使用 hisat2 比对
hisat2 -x genome.hisat2 \
	-1 clean_data/Cs1_rep1_1.fastp.fastq.gz \
	-2 clean_data/Cs1_rep1_2.fastp.fastq.gz \
	--new-summary --rna-strandness RF \
	-S aligned/Cs1_rep1.sam 2> aligned/Cs1_rep1.hisat2.log

# 使用 samtools 排序并构建索引
mamba install -c bioconda samtools
samtools sort -o Cs1_rep1.bam aligned/Cs1_rep1.sam
samtools index Cs1_rep1.bam

# 自动判断链特异性类型
mamba install -c bioconda how_are_we_stranded_here
check_strandedness --gtf ref/genes.gtf \
	--transcripts ref/transcriptome.fa \
	--reads_1 clean_data/Cs2_rep1_1.fastp.fastq.gz \
	--reads_2 clean_data/Cs2_rep1_2.fastp.fastq.gz

# 表达定量
git clone http://git.genek.cn:3333/zhxd2/RunFeatureCounts.git
Rscript software/RunFeatureCounts/run-featurecounts.R \
	-b aligned/Cs1_rep1.bam \
	-g ref/genes.gtf \
	-s 2 \
	-o quanti/Cs1_rep1

# 合并表达矩阵
ls quanti/*.count >quanti_files.txt
perl software/RunFeatureCounts/abundance_estimates_to_matrix.pl \
	--est_method featureCounts \
      	--out_prefix quanti/genes \
      	--quant_files quanti_files.txt

sed -i '1s/^/gene_id/' quanti/genes.counts.matrix
sed -i '1s/^/gene_id/' quanti/genes.TMM.TPM.matrix

# 差异表达分析
perl /pub/software/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/run_DE_analysis.pl \
        --matrix quanti/genes.counts.matrix \
	--method DESeq2 \
        --samples_file samples.txt \
	--contrasts contrasts.txt \
        --output DE_analysis

awk 'FNR==1 && NR!=1{{next}}{{print}}' DE_analysis/*DE_results > quanti/de_result.tsv

sed -i '1s/^/gene_id\t/' quanti/de_result.tsv
rm -rf DE_analysis    

# emapper 功能注释
mamba install -c bioconda -c conda-forge eggnog-mapper
download_eggnog_data.py 
emapper.py -m mmseqs \
	-i ref/proteins.fa \
	-o eggnog \
	--output_dir annotations \
	--cpu 20 \
	--data_dir /pub/database/eggnog-mapper-data

mv annotations/eggnog.emapper.annotations \
	annotations/eggnog_annotations.tsv

# 构建 OrgDB
git clone http://git.genek.cn:3333/zhxd2/emcp.git
Rscript software/emcp/emapperx.R \
	annotations/eggnog_annotations.tsv \
	ref/proteins.fa annotations

