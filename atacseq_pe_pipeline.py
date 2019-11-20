configfile: "config.yml"

# SAMPLES, = glob_wildcards("raw_data/{smp}_R1.fastq.gz")
SAMPLES = [x for condition in config["samples"] for x in config["samples"][condition]]
BT2INDEX = config["bt2_index"]
BLACKLIST = config["blacklist"]
COUNTFILE = [config["countfile"]]
ANNOTATION = config["annotation"]

def revLookup(sample):
    '''Finds the condition grouping name and rep number for a given sample'''
    return next("%s_rep%d" % (key, value.index(sample)) for key, value in config['samples'].items() if sample in value)

mappings = {revLookup(x) : x for x in SAMPLES}

ALL_FASTQC = expand("fastqc_out/{sample}_R1_fastqc.zip", sample=SAMPLES)
ALL_BAMCOV = expand("results/{sample}.dedup.masked.rpkm.bw", sample=list(mappings.keys()))
PEAKS_NARROWPEAK = expand("peaks/{sample}_peaks.narrowPeak", sample=SAMPLES)
ANNOTATED_PEAKS = ["peaks/merged_peaks_annotated.txt"]

rule all:
    input:
        ALL_FASTQC + ALL_BAMCOV + PEAKS_NARROWPEAK + COUNTFILE + ANNOTATED_PEAKS + ["multiqc_report.html"]

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
        "docker://dukegcb/fastqc"
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
        "docker://dukegcb/trim-galore"
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
        "docker://duke-gcb/bowtie2"
    shell:
        "bowtie2 --very-sensitive --no-unal -p {threads} -x {BT2INDEX} -1 {input.READ1} -2 {input.READ2} 2> {log} | samtools view -bS -o {output}"

rule samtools_sort:
    input: rules.bowtie_align.output
        # "results/{sample}.bam"
    output:
        "results/{sample}.sorted.bam"
    log:
        "logs/{sample}.samtools_sort.log"
    singularity:
        "docker://jweinstk/samtools"
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
        "docker://broadinstitute/picard"
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
        "docker://biowardrobe2/bedtools2:v2.26.0"
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
        "docker://jweinstk/samtools"
    shell:
        "samtools index {input} {output} &> {log}"

rule bam_coverage:
    input:
        BAM = lambda wildcards: "results/{}.sorted.dedup.masked.bam".format(mappings[wildcards.sample]),
        BAI = lambda wildcards: "results/{}.sorted.dedup.masked.bai".format(mappings[wildcards.sample])
    output:
        "results/{sample}.dedup.masked.rpkm.bw"
    singularity:
        "docker://genomicpariscentre/deeptools"
    threads: 4
    log:
        "logs/{sample}.bamcoverage.log"
    shell: """
            bamCoverage -p {threads} -b {input.BAM} -o {output} \
                --binSize 10 \
                --normalizeUsingRPKM \
                --ignoreForNormalization chrX chrY chrM\
                --samFlagExclude 1024 \
                &> {log}
                """

rule callpeaks_narrow:
    input:
        "results/{sample}.sorted.dedup.masked.bam"
    output:
        "peaks/{sample}_peaks.narrowPeak"
    singularity:
        "docker://dukegcb/macs2"
    log:
        "logs/{sample}.macs2.log"
    shell:
        "macs2 callpeak --nomodel -t {input} -f BAM \
        --shift -100 --extsize 200 \
        -n peaks/{wildcards.sample} -g mm -q 0.1 &> {log}"

rule idr:
    input:
        lambda wildcards: expand("peaks/{replicate}_peaks.narrowPeak", replicate=config['samples'][wildcards.sample])
    output:
        "peaks/{sample}_idr.txt"
    shell:
        "idr --samples {input} --idr-threshold 0.05 --output-file {output}"

rule combine_pk:
    input:
        [expand("peaks/{sample}_idr.txt", sample=config["samples"]) if config['idr'] else expand("peaks/{sample}_peaks.narrowPeak", sample=SAMPLES)]
    output:
        "peaks/combined_peaks.sorted.bed"
    shell:
        "cat {input} | sort -k1,1 -k2,2n > {output}"

rule merge_peaks:
    input:
        "peaks/combined_peaks.sorted.bed"
    output:
        "peaks/combined_peaks_merged.bed"
    singularity:
        "docker://biowardrobe2/bedtools2:v2.26.0"
    shell:
        """bedtools merge -i {input} > {output}"""

rule bed_to_saf:
    input:
        "peaks/combined_peaks_merged.bed"
    output:
        "peaks/combined_peaks_merged.saf"
    shell:
        """
        awk 'OFS="\\t" {{print $1"."$2+1"."$3, $1, $2+1, $3, "."}}' {input} > {output}
        """

rule featureCounts:
    input:
        FILES=expand("results/{sample}.sorted.dedup.masked.bam", sample=SAMPLES),
        PEAKSET="peaks/combined_peaks_merged.saf"
    output:
        COUNTFILE
    threads: 4
    singularity:
        "docker://genomicpariscentre/featurecounts"
    log:
        "logs/featureCounts.log"
    shell:
        "featureCounts -T {threads} --largestOverlap -F SAF --ignoreDup \
        -a {input.PEAKSET} \
        -o {output} \
        {input.FILES} \
        &> {log}"

rule saf_to_bed:
    input:
        "peaks/combined_peaks_merged.saf"
    output:
        temp("peaks/combined_peaks_merged_named.bed")
    shell:
        """
        awk 'OFS="\\t" {{print $2,$3-1,$4,$1,$5}}' {input} > {output}
        """

rule annotate_peaks:
    input:
        "peaks/combined_peaks_merged_named.bed"
    output:
        "peaks/merged_peaks_annotated.txt"
    shell:
        "annotatePeaks.pl {input} mm10 -gtf {ANNOTATION} > {output}"

rule multiqc:
    input:
        ALL_FASTQC + ALL_BAMCOV + PEAKS_NARROWPEAK + COUNTFILE + ANNOTATED_PEAKS
    output: "multiqc_report.html"
    singularity:
        "docker://ewels/multiqc"
    threads: 1
    shell:
        "multiqc ."
