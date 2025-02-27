import pandas as pd

# 指定配置文件  
configfile: "config.yaml"

# 通过配置文件读取样本列表
SAMPLES = pd.read_csv(config["samples_file"], sep='\t', header=None, usecols=[1]).squeeze("columns").tolist()

rule all:
  input:
    expand("aligned/{sample}.bam", sample=SAMPLES),
    "quanti/genes.counts.matrix",
    "quanti/genes.TMM.TPM.matrix",
    "quanti/de_result.tsv",
    "clean_data/fastp_summary.tsv",
    "aligned/hisat2_summary.tsv",
    "annotations/eggnog_annotations.tsv",
    "annotations/org.My.eg.db_1.0.tar.gz"

rule trim_fastq:
  input:
    r1="raw_data/{sample}_1.fastq.gz",
    r2="raw_data/{sample}_2.fastq.gz"
  output:
    r1="clean_data/{sample}_1.fastp.fastq.gz",
    r2="clean_data/{sample}_2.fastp.fastq.gz",
    report="clean_data/{sample}.fastp.html",
    json="clean_data/{sample}.fastp.json"
  threads: 4
  log: "logs/{sample}.fastp.log"
  conda: "envs/fastp.yaml"
  shell:
    """
      fastp -i {input.r1} -I {input.r2} \
        -o {output.r1} -O {output.r2} \
        -h {output.report} -j {output.json} \
        -w {threads} 1> {log} 2>&1
    """

rule combine_fastp_report:
  input:
    expand("clean_data/{sample}.fastp.json", sample=SAMPLES)
  output:
    "clean_data/fastp_summary.tsv"
  shell:
    """
      Rscript scripts/combine_fastp_report.R -i clean_data/ -o clean_data/
    """

rule build_hisat2_index:
  input:
    genome=config["genome"]
  output:
    idx=expand("genome.hisat2.{idx}.ht2", idx=range(1, 9))
  log: "logs/build_hisat2_index.log"
  conda: "envs/hisat2.yaml"
  shell:
    """
      hisat2-build {input.genome} genome.hisat2 1>{log} 2>&1
    """

rule align_hisat2:
  input:
    r1="clean_data/{sample}_1.fastp.fastq.gz",
    r2="clean_data/{sample}_2.fastp.fastq.gz",
    idx=expand("genome.hisat2.{idx}.ht2", idx=range(1, 9))
  output:
    sam=temp("aligned/{sample}.sam"),
    log="aligned/{sample}.hisat2.log"
  params:
    libtype=config["libtype"]
  threads: 4
  conda: "envs/hisat2.yaml"
  shell:
    """
      hisat2 -x genome.hisat2 -1 {input.r1} -2 {input.r2} \
        --new-summary --rna-strandness RF \
        -S {output.sam} \
        1> {output.log} 2>&1
    """

rule combine_hisat2_report:
  input:
    expand("aligned/{sample}.hisat2.log", sample=SAMPLES)
  output:
    "aligned/hisat2_summary.tsv"
  shell:
    """
      python scripts/combine_hisat2_report.py -i aligned/ -o aligned/
    """

rule samtools_sort_idx:
  input:
    sam="aligned/{sample}.sam"
  output:
    bam=protected("aligned/{sample}.bam"),
    idx="aligned/{sample}.bam.bai"
  threads: 4
  log: "logs/{sample}.sort.log"
  conda: "envs/samtools.yaml"
  shell:
    """
      samtools sort -@ {threads} -o {output.bam} {input.sam} 1>{log} 2>&1
      samtools index {output.bam}
    """

hisat2_strandness = config["libtype"]
def convert_strandness(hisat2_strandness):
    if hisat2_strandness == "FR":
        return 1
    elif hisat2_strandness == "RF":
        return 2
    else:
        raise ValueError("Invalid HISAT2 strandness value")

rule quanti_featureCounts:
  input:
    bam="aligned/{sample}.bam",
    gtf="ref/genes.gtf"
  output:
    expr="quanti/{sample}.count",
    log="quanti/{sample}.log"
  conda: "envs/RunFeatureCounts.yaml"
  params:
    strandSpecific=convert_strandness(config["libtype"]),
    prefix="quanti/{sample}"
  shell:
    """
      Rscript software/RunFeatureCounts/run-featurecounts.R \
	-b {input.bam} -g {input.gtf} -s {params.strandSpecific} \
	-o {params.prefix}
    """

rule abundance_estimates_to_matrix:
  input:
    exprs=expand("quanti/{sample}.count", sample=SAMPLES)
  output:
    counts_mat="quanti/genes.counts.matrix",
    tpm_mat="quanti/genes.TMM.TPM.matrix"
  log: "logs/abundance_estimates_to_matrix.log"
  conda: "envs/trinity.yaml"
  shell:
    """
    ls {input.exprs} >quant_files.txt
    perl software/RunFeatureCounts/abundance_estimates_to_matrix.pl \
	--est_method featureCounts \
        --out_prefix quanti/genes \
        --quant_files quant_files.txt 1>{log} 2>&1

    # add feature id
    sed -i '1s/^/gene_id/' {output.counts_mat}
    sed -i '1s/^/gene_id/' {output.tpm_mat}
    """

rule de_analysis:
  input:
    counts_mat="quanti/genes.counts.matrix"
  output:
    de_result="quanti/de_result.tsv"
  params:
    TRINITY_HOME=config["TRINITY_HOME"],
    method=config["DE_method"],
    samples_file=config["samples_file"],
    contrasts_file=config["contrasts_file"]
  log: "logs/de_analysis.log"
  conda: "envs/trinity.yaml"
  shell:
    """
      perl {params.TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl \
        --matrix {input.counts_mat} --method {params.method} \
        --samples_file {params.samples_file} --contrasts {params.contrasts_file} \
        --output DE_analysis

    # Combine all DE_results files
      awk 'FNR==1 && NR!=1{{next}}{{print}}' DE_analysis/*DE_results > {output.de_result}

    # Remove the diff_expr_results directory
      rm -rf DE_analysis
    """

rule emapper:
  input:
    protein=config["protein"]
  output:
    result="annotations/eggnog_annotations.tsv"
  threads: 40
  params:
    emapper_mode=config["emapper_mode"],
    eggnog_db=config["eggnog_db"]
  shell:
    """
    emapper.py -m {params.emapper_mode} \
        -i {input.protein} \
        -o eggnog \
        --output_dir annotations \
        --cpu {threads} \
        --data_dir {params.eggnog_db}

    mv annotations/eggnog.emapper.annotations {output.result}
    """

rule create_orgdb:
  input:
    emapper_anno="annotations/eggnog_annotations.tsv",
    protein=config["protein"]
  output:
    "annotations/org.My.eg.db_1.0.tar.gz"
  shell:
    """
    Rscript software/emcp/emapperx.R {input.emapper_anno} \
        {input.protein} annotations
    """
