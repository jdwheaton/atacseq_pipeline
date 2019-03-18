configfile: "config.yml"

SAMPLES, = glob_wildcards("raw_data/{smp}_R1.fastq.gz")
BT2INDEX = config["bt2_index"]
BLACKLIST = config["blacklist"]

ALL_FASTQC = expand("fastqc_out/{sample}_R1_fastqc.zip", sample=SAMPLES)
ALL_BAMCOV = expand("results/{sample}.dedup.masked.rpkm.bw", sample=SAMPLES)

rule all:
    input:
        ALL_FASTQC + ALL_BAMCOV

rule fastqc:
    input:
        "raw_data/{sample}_R1.fastq.gz",
        "raw_data/{sample}_R2.fastq.gz"
    output:
        "fastqc_out/{sample}_R1_fastqc.html",
        "fastqc_out/{sample}_R1_fastqc.zip",
        "fastqc_out/{sample}_R2_fastqc.html",
        "fastqc_out/{sample}_R2_fastqc.zip"
    log:
        "logs/{sample}.fastqc.log"
    singularity:
        "shub://jdwheaton/singularity-ngs:qc_trim"
    threads: 1
    shell:
        "fastqc -o fastqc_out/ {input} &> {log}"

rule trim_galore:
    input:
        "raw_data/{sample}_R1.fastq.gz",
        "raw_data/{sample}_R2.fastq.gz"
    output:
        "trimmed_fastq/{sample}_R1_val_1.fq",
        "trimmed_fastq/{sample}_R2_val_2.fq"
    log:
        "logs/{sample}.trimgalore.log"
    threads: 1
    singularity:
        "shub://jdwheaton/singularity-ngs:qc_trim"
    shell:
        "trim_galore --dont_gzip --fastqc -o trimmed_fastq/ --paired {input} &> {log}"

rule bowtie_align:
    output:
        "results/{sample}.bam"
    input:
        READ1 = "trimmed_fastq/{sample}_R1_val_1.fq",
        READ2 = "trimmed_fastq/{sample}_R2_val_2.fq"
    log:
        "logs/{sample}.bowtie.log"
    threads: 5
    singularity:
        "shub://jdwheaton/singularity-ngs:latest"
    shell:
        "bowtie2 --very-sensitive --no-unal -p {threads} -x {BT2INDEX} -1 {input.READ1} -2 {input.READ2} 2> {log} | samtools view -bS -o {output}"

rule samtools_sort:
    input:
        "results/{sample}.bam"
    output:
        "results/{sample}.sorted.bam"
    log:
        "logs/{sample}.samtools_sort.log"
    singularity:
        "shub://jdwheaton/singularity-ngs:latest"
    shell:
        "samtools sort -o {output} {input} &> {log}"

rule remove_duplicates:
    input:
        "results/{sample}.sorted.bam"
    output:
        BAM = "results/{sample}.sorted.dedup.bam",
        METRICS = "results/{sample}.dedup.metrics.txt"
    log:
        "logs/{sample}.rmdup.log"
    singularity:
        "shub://jdwheaton/singularity-ngs:chip_atac_post"
    shell: """
            java -Xmx16g -jar /picard.jar MarkDuplicates \
                INPUT={input} \
                OUTPUT={output.BAM} \
                METRICS_FILE={output.METRICS} \
                ASSUME_SORT_ORDER=coordinate \
                REMOVE_DUPLICATES=false \
                2> {log}
                """

rule remove_blacklisted:
    input:
        BAM = "results/{sample}.sorted.dedup.bam",
        BLACKLIST = BLACKLIST
    output:
        "results/{sample}.sorted.dedup.masked.bam"
    log:
        "logs/{sample}.deblacklist.log"
    singularity:
        "shub://jdwheaton/singularity-ngs:chip_atac_post"
    shell:
        "bedtools intersect -v -a {input.BAM} -b {input.BLACKLIST} > {output} 2> {log}"

rule samtools_index:
    input:
        "results/{sample}.sorted.dedup.masked.bam"
    output:
        "results/{sample}.sorted.dedup.masked.bai"
    log:
        "logs/{sample}.sami_ndex.log"
    singularity:
        "shub://jdwheaton/singularity-ngs:latest"
    shell:
        "samtools index {input} {output} &> {log}"

rule bam_coverage:
    input:
        BAM = "results/{sample}.sorted.dedup.masked.bam",
        BAI = "results/{sample}.sorted.dedup.masked.bai"
    output:
        "results/{sample}.dedup.masked.rpkm.bw"
    singularity:
        "shub://jdwheaton/singularity-ngs:chip_atac_post"
    threads: 4
    log:
        "logs/{sample}.bamcoverage.log"
    shell: """
            bamCoverage -p {threads} -b {input.BAM} -o {output} \
                --binSize 10 \
                --normalizeUsing RPKM \
                --ignoreForNormalization chrX chrY chrM\
                --samFlagExclude 1024 \
                &> {log}
                """
