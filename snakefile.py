import pandas as pd
import os
import glob
import sys

configfile: "config.yaml"

# ------------------------------------------
# Load configuration parameters
# ------------------------------------------
DATA_DIR       = config["data_dir"]
FASTQ_PATTERN  = config["fastq_pattern"]
GENOMES        = config["genomes"]
QVALS          = config["qvals"]
READS          = config["reads"]
SAMPLES_FILTER = config["samples"]
RESOURCES_DIR  = config["resources_dir"]
LOG_DIR = config["log_dir"]
RESULTS_DIR = config["results_dir"]
MAX_ALIGNMENTS = config["k"]
# ------------------------------------------
# Discover samples automatically
# ------------------------------------------
samples_raw, reads_raw = glob_wildcards(os.path.join(DATA_DIR, FASTQ_PATTERN))
SAMPLES = sorted(set(samples_raw))
if SAMPLES_FILTER:
    SAMPLES = [s for s in SAMPLES if s in SAMPLES_FILTER]

print(f"→ DATA_DIR:       {DATA_DIR}")
print(f"→ RESOURCES_DIR:  {RESOURCES_DIR}")
print(f"→ SAMPLES:        {SAMPLES}")
print(f"→ READS:          {READS}")
print(f"→ MAX ALIGNMENTS: {MAX_ALIGNMENTS}")

# ------------------------------------------
# Rule: all
# ------------------------------------------
rule all:
    input:
        # Genome FASTA and chromosome sizes files
        expand(
            f"{RESOURCES_DIR}/{{genome}}{{suffix}}",
            genome=GENOMES, suffix=['.fa', '.sizes.txt']
        ),
        
        # Genome indexes
        expand(
            f"{RESOURCES_DIR}/{{genome}}.fa{{suffix}}",
            genome=GENOMES,
            suffix=[
                ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                ".rev.1.bt2", ".rev.2.bt2"
            ]
        ),

        # FastQC (raw)
        expand(
            f"{RESULTS_DIR}/qc/fastqc_raw/{{sample}}_{{read}}_fastqc.html",
            sample=SAMPLES,
            read=READS
        ),

        # Cutadapt outputs
        expand(
            f"{RESULTS_DIR}/cutadapt/{{sample}}_R1.fastq.gz",
            sample=SAMPLES
        ),

        # FastQC (trimmed)
        expand(
            f"{RESULTS_DIR}/qc/fastqc_cutadapt/{{sample}}_{{read}}_fastqc.html",
            sample=SAMPLES,
            read=READS
        ),

        # Alignment
        expand(
            f"{RESULTS_DIR}/bowtie2_align/{{sample}}-{{genome}}-{{k}}.bam",
            sample=SAMPLES,
            genome=GENOMES, 
            k=MAX_ALIGNMENTS
        ),

        # Filtering
        expand(
            f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}-{{k}}.bam",
            sample=SAMPLES,
            genome=GENOMES,
            k=MAX_ALIGNMENTS
        ),

        # MarkDuplicates
        expand(
            f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}-{{k}}.bam",
            sample=SAMPLES,
            genome=GENOMES,
            k=MAX_ALIGNMENTS
        ),

        # BEDPE & tagAlign
        expand(
            f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}-{{k}}.bedpe",
            sample=SAMPLES,
            genome=GENOMES,
            k=MAX_ALIGNMENTS
        ),

        # MACS3 signal bedGraph
        expand(
            f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}-{{k}}/{{sample}}-{{genome}}-{{k}}_treat_pileup.bdg",
            sample=SAMPLES,
            genome=GENOMES,
            k=MAX_ALIGNMENTS
        ),

        # MACS3 callpeak (but only to obtain signal bigWig)
        expand(
            f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}-{{k}}/{{sample}}-{{genome}}-{{k}}.bigWig",
            sample=SAMPLES,
            genome=GENOMES,
            k=MAX_ALIGNMENTS
        ),

        # MACS3 callpeak
        expand(
            f"{RESULTS_DIR}/macs3_callpeak/{{sample}}-{{genome}}-{{k}}-{{q}}/{{sample}}-{{genome}}-{{k}}-{{q}}_peaks.narrowPeak",
            sample=SAMPLES,
            genome=GENOMES,
            k=MAX_ALIGNMENTS,
            q=QVALS
        ),

        # Fragment length distribution
        expand(
            f"{RESULTS_DIR}/qc/frag_len/{{sample}}-{{genome}}-{{k}}.png",
            sample=SAMPLES,
            genome=GENOMES,
            k=MAX_ALIGNMENTS
        ),

        # Fragment length distribution
        expand(
            f"{RESULTS_DIR}/qc/tss_enrichment/{{sample}}-{{genome}}-{{k}}.png",
            sample=SAMPLES,
            genome=GENOMES,
            k=MAX_ALIGNMENTS
        )
        # # BigWig
        # expand(
        #     "results/macs3/{sample}-{genome}/{sample}-{genome}.bigWig",
        #     sample=SAMPLES,
        #     genome=GENOMES
        # ),

        # # QC: fragment length
        # expand(
        #     "results/qc/frag_len/{sample}-{genome}.png",
        #     sample=SAMPLES,
        #     genome=GENOMES
        # ),
        # expand(
        #     "results/qc/frag_len/{sample}-{genome}.txt",
        #     sample=SAMPLES,
        #     genome=GENOMES
        # ),

        # # QC: TSS enrichment
        # expand(
        #     "results/qc/tss_enrichment/{sample}-{genome}.png",
        #     sample=SAMPLES,
        #     genome=GENOMES
        # ),
        # expand(
        #     "results/qc/tss_enrichment/{sample}-{genome}.txt",
        #     sample=SAMPLES,
        #     genome=GENOMES
        # ),

        # # QC: FRiP
        # expand(
        #     "results/qc/frip_all/{sample}-{genome}-{q}.txt",
        #     sample=SAMPLES,
        #     genome=GENOMES,
        #     q=QVALS
        # ),

        # # ChromBPNet prep & bias
        # expand(
        #     "results/chrombpnet/{sample}-{genome}/chrombpnet_prep_nonpeaks.ok",
        #     sample=SAMPLES,
        #     genome=GENOMES
        # ),
        # expand(
        #     "results/chrombpnet/{sample}-{genome}/chrombpnet_bias.ok",
        #     sample=SAMPLES,
        #     genome=GENOMES
        # ),

        # # Pseudoreplicates & IDR
        # expand(
        #     "results/pseudoreps/{sample}.pseudorep1.tagalign",
        #     sample=SAMPLES
        # ),
        # expand(
        #     "results/pseudoreps/{sample}.pseudorep2.tagalign",
        #     sample=SAMPLES
        # ),
        # expand(
        #     "results/macs3_pseudoreps/{sample}-{pseudorep}/{sample}-{pseudorep}_peaks.narrowPeak",
        #     sample=SAMPLES,
        #     pseudorep=["pseudorep1","pseudorep2"]
        # ),
        # expand(
        #     "results/idr/{sample}.IDR.bed",
        #     sample=SAMPLES
        # ),
        # expand(
        #     "results/overlap_peaks/{sample}.overlap_peaks.bed",
        #     sample=SAMPLES
        # )

# ------------------------------------------
# Example: get_genome rule using RESOURCES_DIR
# ------------------------------------------
rule get_genome:
    output:
        fa    = f"{RESOURCES_DIR}/{{genome}}.fa",
        sizes = f"{RESOURCES_DIR}/{{genome}}.sizes.txt"
    threads: 1
    log: "logs/get_genome/{genome}.log"
    conda: "MOHD-ATAC"
    resources:
            slurm_partition='4hours', runtime=30, mem_mb='4G', cpus_per_task=1

    shell:
        """
        (
        wget -q -O {output.fa}.gz \
            $( \
              if [ {wildcards.genome} = 'GRCh38' ]; then \
                echo 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'; \
              else \
                echo 'http://t2t.gi.ucsc.edu/chm13/.../t2t-chm13-v1.0.fa.gz'; \
              fi \
            ) \
        && gzip -d -c {output.fa}.gz > {output.fa} \
        && rm {output.fa}.gz \
        && samtools faidx {output.fa} \
        && cut -f1,2 {output.fa}.fai > {output.sizes}
        ) &> {log}
        """

rule bowtie2_build:
    input:
        fa = f"{RESOURCES_DIR}/{{genome}}.fa"
    output:
        multiext(
            f"{RESOURCES_DIR}/{{genome}}.fa",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
            ".rev.1.bt2", ".rev.2.bt2"
        )
    resources:
        slurm_partition='4hours', runtime=240, mem_mb='96G', cpus_per_task=24

    threads: 24
    log: f"{LOG_DIR}/bowtie2_build/{{genome}}.log"
    conda: "MOHD-ATAC"
    shell:
        """
        (
            echo "Building Bowtie2 index for {wildcards.genome}"
            bowtie2-build --threads {threads} {input.fa} {RESOURCES_DIR}/{wildcards.genome}.fa
        ) &> {log}
        """


    
rule fastqc_raw:
    input: 
        os.path.join(DATA_DIR, "{sample}_{read}.fastq.gz")
    output: 
        f"{RESULTS_DIR}/qc/fastqc_raw/{{sample}}_{{read}}_fastqc.html"
    conda:
        "MOHD-ATAC"
    threads: 8
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: f"{LOG_DIR}/fastqc_raw/{{sample}}-{{read}}.log"
    params: 
        o = os.path.join(RESULTS_DIR, 'qc', 'fastqc_raw')
    shell:
        """
        (
        fastqc -t {threads} -o {params.o} {input}
        ) &> {log}
        """

rule cutadapt:
    input: 
        R1 = os.path.join(DATA_DIR, "{sample}_R1.fastq.gz"),
        R2 = os.path.join(DATA_DIR, "{sample}_R2.fastq.gz")
    output:
        R1 = f"{RESULTS_DIR}/cutadapt/{{sample}}_R1.fastq.gz",
        R2 = f"{RESULTS_DIR}/cutadapt/{{sample}}_R2.fastq.gz"
    threads: 8
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    conda:
        "MOHD-ATAC"
    log:
        f"{LOG_DIR}/cutadapt/{{sample}}.log",
    
    shell:
        """
        (
        cutadapt -q 10 -m 15 -e 0.10 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -j {threads} -o {output.R1} -p {output.R2} {input.R1} {input.R2}
        ) &> {log}
        """
        
        
rule fastqc_cutadapt:
    input: 
        f"{RESULTS_DIR}/cutadapt/{{sample}}_{{read}}.fastq.gz"
    output: 
        f"{RESULTS_DIR}/qc/fastqc_cutadapt/{{sample}}_{{read}}_fastqc.html"
    conda:
        "MOHD-ATAC"
    threads: 8
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: 
        f"{LOG_DIR}/fastqc_cutadapt/{{sample}}-{{read}}.log"
    params: 
        o = os.path.join(RESULTS_DIR, 'qc', 'fastqc_cutadapt')
    shell:
        """
        (
        fastqc -t {threads} -o {params.o} {input}
        ) &> {log}
        """
        
rule bowtie2_align:
    input:
        fa = f"{RESOURCES_DIR}/{{genome}}.fa",

        index = multiext(
            f"{RESOURCES_DIR}/{{genome}}.fa",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
            ".rev.1.bt2", ".rev.2.bt2"
        ),

        R1 = f"{RESULTS_DIR}/cutadapt/{{sample}}_R1.fastq.gz",
        R2 = f"{RESULTS_DIR}/cutadapt/{{sample}}_R2.fastq.gz"

    output:
        f"{RESULTS_DIR}/bowtie2_align/{{sample}}-{{genome}}-{{k}}.bam"
    threads: 8 
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log:
        f"{LOG_DIR}/bowtie2_align/{{sample}}-{{genome}}-{{k}}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        echo "Running bowtie2 align"
        bowtie2 -X 2000 --mm -k {wildcards.k} \
        --threads {threads} \
        -x {input.fa} \
        --rg-id {wildcards.sample}-{wildcards.genome}-{wildcards.k} \
        --rg SM:{wildcards.sample}-{wildcards.genome}-{wildcards.k} \
        -1 {input.R1} \
        -2 {input.R2} | \
        samtools view -1 -S /dev/stdin > /tmp/{wildcards.sample}-{wildcards.genome}-{wildcards.k}.tmp.bam
        samtools sort -@ {threads} -T {wildcards.sample}-{wildcards.genome}-{wildcards.k} -o {output} /tmp/{wildcards.sample}-{wildcards.genome}-{wildcards.k}.tmp.bam
        samtools flagstat -@ {threads} {output} > {output}.flagstat
        rm /tmp/{wildcards.sample}-{wildcards.genome}-{wildcards.k}.tmp.bam
        ) &> {log}
        """

def get_exclude_flag(wildcards):
    if int(wildcards.k) > 1:
        return("1548")
    else:
        return("1804")
        
rule filter:
    input: 
        f"{RESULTS_DIR}/bowtie2_align/{{sample}}-{{genome}}-{{k}}.bam"
    output: 
        f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}-{{k}}.bam"
    threads: 8 
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: 
        f"{LOG_DIR}/filter/{{sample}}-{{genome}}-{{k}}.log"
    conda:
        "MOHD-ATAC"
    params:
        exclude_flag = lambda wildcards: get_exclude_flag(wildcards),
        prefix = lambda wildcards: "-".join(wildcards)
    shell:
        """
        (
        echo {params.prefix}
        samtools view -F {params.exclude_flag} -f 2 -q 30 -u {input} | \
        samtools sort -n /dev/stdin -o /tmp/{params.prefix}.tmp.nmsrt.bam -T {params.prefix} -@ {threads}

        samtools view -h /tmp/{params.prefix}.tmp.nmsrt.bam |
        samtools fixmate -r /dev/stdin /tmp/{params.prefix}.tmp.fixmate.bam
        
        samtools view -F {params.exclude_flag} -f 2 -u /tmp/{params.prefix}.tmp.fixmate.bam | \
            samtools sort /dev/stdin -@ {threads} -o /tmp/{params.prefix}.tmp.coordsort.bam 

        samtools view -h /tmp/{params.prefix}.tmp.coordsort.bam | grep -v "chrM" | samtools view -b > {output}
        
        samtools index -@ {threads} {output} {output}.bai
        samtools flagstat -@ {threads} {output} > {output}.flagstat

        rm /tmp/{params.prefix}.tmp.nmsrt.bam
        rm /tmp/{params.prefix}.tmp.fixmate.bam
        rm /tmp/{params.prefix}.tmp.coordsort.bam
        ) &> {log}
        """

rule picard:
    input: 
        f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}-{{k}}.bam"
    output: 
        f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}-{{k}}.bam", 
        f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}-{{k}}.txt"
    threads: 8 
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: 
        f"{LOG_DIR}/picard/{{sample}}-{{genome}}-{{k}}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        picard MarkDuplicates \
            INPUT={input} \
            OUTPUT={output[0]} \
            METRICS_FILE={output[1]} \
            REMOVE_DUPLICATES=TRUE \
            ASSUME_SORT_ORDER=coordinate \
            USE_JDK_DEFLATER=TRUE \
            USE_JDK_INFLATER=TRUE \
            VALIDATION_STRINGENCY=LENIENT
        samtools index -@ {threads} {output[0]} {output[0]}.bai
        samtools flagstat -@ {threads} {output[0]} > {output[0]}.flagstat
        ) &> {log}
        """

rule bam_to_tagalign:
    input: 
        f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}-{{k}}.bam"
    output: 
        bedpe = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}-{{k}}.bedpe",
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}-{{k}}.tagalign.gz"
            
    threads: 8 
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: 
        f"{LOG_DIR}/bedpe_tagalign/{{sample}}-{{genome}}-{{k}}.log"
    conda:
        "MOHD-ATAC"
    params:
            prefix = lambda wildcards: "-".join(wildcards),
            results_dir = os.path.join(RESULTS_DIR, 'bedpe_tagalign')
    shell:
        """
        (
        python scripts/bam_to_tagalign.py {input} {output.bedpe} {params.results_dir}/{params.prefix}.tagalign --threads {threads}
        touch {output.tagalign}
        ) &> {log}
        """

rule macs3_signal:
    input:
        f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}-{{k}}.tagalign.gz"
    output:
        f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}-{{k}}/{{sample}}-{{genome}}-{{k}}_treat_pileup.bdg"
    threads: 8 
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    conda:
        "MOHD-ATAC"
    log: 
        f"{LOG_DIR}/macs3_signal/{{sample}}-{{genome}}-{{k}}.log"
    params:
            prefix = lambda wildcards: "-".join(wildcards),
            results_dir = os.path.join(RESULTS_DIR, 'macs3_signal')
    shell:
        """
        (
        macs3 callpeak -f BED -t {input} -n {params.prefix} --outdir {params.results_dir}/{params.prefix} --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q .01 -B
        touch {output}
        ) &> {log}
        """

rule bedGraphToBigWig:
    input: 
        f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}-{{k}}/{{sample}}-{{genome}}-{{k}}_treat_pileup.bdg"
    output:
        f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}-{{k}}/{{sample}}-{{genome}}-{{k}}.bigWig"
    threads: 1
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='4G', cpus_per_task=1
    conda:
        "MOHD-ATAC"
    log:
        f"{LOG_DIR}/bedGraphToBigWig/{{sample}}-{{genome}}-{{k}}.log"
    params:
            prefix = lambda wildcards: "-".join(wildcards),
            results_dir = os.path.join(RESULTS_DIR, 'macs3_signal')
    shell:
        """
        (
        echo "Fixing bedGraph"
        bedtools intersect -a {input} -b resources/{wildcards.genome}.sort.bed -sorted > {params.results_dir}/{params.prefix}/{params.prefix}.bg
        echo "Running bedGraphToBigWig"
        bedGraphToBigWig {params.results_dir}/{params.prefix}/{params.prefix}.bg resources/{wildcards.genome}.sizes.txt {output}
        ) &> {log}
        """


rule macs3_callpeak:
    input:
        f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}-{{k}}.tagalign.gz"
    output:
        f"{RESULTS_DIR}/macs3_callpeak/{{sample}}-{{genome}}-{{k}}-{{q}}/{{sample}}-{{genome}}-{{k}}-{{q}}_peaks.narrowPeak"
    threads: 8 
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    conda:
        "MOHD-ATAC"
    log: 
        f"{LOG_DIR}/macs3_callpeak/{{sample}}-{{genome}}-{{k}}-{{q}}.log"
    params:
            prefix = lambda wildcards: "-".join(wildcards),
            results_dir = os.path.join(RESULTS_DIR, 'macs3_callpeak')

    shell:
        """
        (
        macs3 callpeak -f BED -t {input} -n {params.prefix} --outdir {params.prefix}/{params.prefix} --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q {wildcards.q}
        touch {output}
        ) &> {log}
        """

rule frag_len:
    input:
        f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}-{{k}}.bam"
    output: 
        f"{RESULTS_DIR}/qc/frag_len/{{sample}}-{{genome}}-{{k}}.png", 
        f"{RESULTS_DIR}/qc/frag_len/{{sample}}-{{genome}}-{{k}}.txt"
    threads: 8
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: 
        f"{LOG_DIR}/frag_len/{{sample}}-{{genome}}-{{k}}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        python scripts/plot_fragment_length_distr.py {input} {threads} {output[0]} {output[1]}
        ) &> {log}
        """


rule tss_enrichment:
    input: 
        f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}-{{k}}.bam"
    output:
        f"{RESULTS_DIR}/qc/tss_enrichment/{{sample}}-{{genome}}-{{k}}.png", 
        f"{RESULTS_DIR}/qc/tss_enrichment/{{sample}}-{{genome}}-{{k}}.txt"
    threads: 8 
    resources:
            slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: 
        f"{LOG_DIR}/frag_len/{{sample}}-{{genome}}-{{k}}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        python scripts/tss_enrichment.py {input} resources/ENCFF493CCB.bed.gz {output[0]} {output[1]}
        ) &> {log}
        """

rule frip_all:
    input: 
        peaks = f"{RESULTS_DIR}/macs3_callpeak/{{sample}}-{{genome}}-{{k}}-{{q}}/{{sample}}-{{genome}}-{{k}}-{{q}}_peaks.narrowPeak", 
        reads = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}-{{k}}.tagalign.gz"
    output:
        f"{RESULTS_DIR}/qc/frip/{{sample}}-{{genome}}-{{k}}.txt"
    threads: 8 
    resources:
        slurm_partition='4hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: 
        f"{LOG_DIR}/frip/{{sample}}-{{genome}}-{{k}}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        python scripts/calculate_frip.py {input.reads} {input.peaks} {output}
        ) &> {log}
        """ 


# rule pseudoreps:
#     input: "results/bedpe_tagalign/{sample}.bedpe"
#     output: "results/pseudoreps/{sample}.pseudorep1.tagalign", "results/pseudoreps/{sample}.pseudorep2.tagalign"
#     threads: 1
#     conda:
#         "MOHD-ATAC"
#     log: "logs/pseudoreps/{sample}.log"
#     shell:
#         """
#         (
#         python scripts/split_bedpe_to_replicates.py {input} {output[0]} {output[1]}
#         ) &> {log}
#         """    

# rule macs3_pseudoreps:
#     input: "results/pseudoreps/{sample}.{pseudorep}.tagalign"
#     output:
#         "results/macs3_pseudoreps/{sample}-{pseudorep}/{sample}-{pseudorep}_peaks.narrowPeak"
#     threads: 24
#     conda:
#         "MOHD-ATAC"
#     log:"logs/macs3_pseudoreps/{sample}-{pseudorep}.log"
#     shell:
#         """
#         (
#         macs3 callpeak -f BED -t {input} -n {wildcards.sample}-{wildcards.pseudorep} --outdir results/macs3_pseudoreps/{wildcards.sample}-{wildcards.pseudorep} --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q .01 -B
#         touch {output}
#         ) &> {log}
#         """


# rule idr:
#     input: "results/macs3_pseudoreps/{sample}-pseudorep1/{sample}-pseudorep1_peaks.narrowPeak", "results/macs3_pseudoreps/{sample}-pseudorep2/{sample}-pseudorep2_peaks.narrowPeak"
#     output: "results/idr/{sample}.IDR.bed", "results/idr/{sample}.IDR.png"
#     threads: 4
#     conda:
#         "idr"
#     log: "logs/idr/{sample}.log"
#     shell:
#         """
#         (
#         sort -k8,8nr {input[0]} > {input[0]}.sorted
#         sort -k8,8nr {input[1]} > {input[1]}.sorted
#         idr --samples {input[0]}.sorted {input[1]}.sorted --input-file-type narrowPeak --rank p.value --output-file /tmp/{wildcards.sample}.IDR.unsorted.bed --plot
#         cat /tmp/{wildcards.sample}.IDR.unsorted.bed | sort -k1,1 -k2,2n > {output[0]}
#         mv /tmp/{wildcards.sample}.IDR.unsorted.bed.png {output[1]}
#         ) &> {log}
#         """

# rule overlap_peaks:
#     input: "results/macs3_pseudoreps/{sample}-pseudorep1/{sample}-pseudorep1_peaks.narrowPeak", 
#             "results/macs3_pseudoreps/{sample}-pseudorep2/{sample}-pseudorep2_peaks.narrowPeak"
#     output: "results/overlap_peaks/{sample}.overlap_peaks.bed"
#     threads: 4
#     conda:
#         "MOHD-ATAC"
#     log: "logs/overlap_peaks/{sample}.log"
#     shell:
#         """
#         (
#         bedtools intersect -a {input[0]} -b {input[1]} -u | sort -k1,1 -k2,2n > {output}
#         ) &> {log}
#         """

# rule frag_len:
#     input: "results/picard/{sample}-{genome}.bam"
#     output: "results/qc/frag_len/{sample}-{genome}.png", "results/qc/frag_len/{sample}-{genome}.txt"
#     threads: 8
#     log: "logs/qc/frag_len/{sample}-{genome}.log"
#     conda:
#         "MOHD-ATAC"
#     shell:
#         """
#         (
#         python scripts/plot_fragment_length_distr.py {input} {threads} {output[0]} {output[1]}
#         ) &> {log}
#         """


# rule tss_enrichment:
#     input: "results/picard/{sample}-{genome}.bam"
#     output: "results/qc/tss_enrichment/{sample}-{genome}.png", "results/qc/tss_enrichment/{sample}-{genome}.txt"
#     threads: 8 
#     resources:
#             slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
#     log: "logs/qc/tss_enrichment/{sample}-{genome}.log"
#     conda:
#         "MOHD-ATAC"
#     shell:
#         """
#         (
#         python scripts/tss_enrichment.py {input} resources/ENCFF493CCB.bed.gz {output[0]} {output[1]}
#         ) &> {log}
#         """

# rule frip_all:
#     input: peaks = "results/macs3_callpeak/{sample}-{genome}-{q}/{sample}-{genome}-{q}_peaks.narrowPeak", reads = "results/bedpe_tagalign/{sample}-{genome}.tagalign.gz"
#     output: "results/qc/frip_all/{sample}-{genome}-{q}.txt"
#     threads: 8 
#     resources:
#         slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
#     log: "logs/qc/frip_all/{sample}-{genome}-{q}.log"
#     conda:
#         "MOHD-ATAC"
#     shell:
#         """
#         (
#         python scripts/calculate_frip.py {input.reads} {input.peaks} {output}
#         ) &> {log}
#         """ 

# rule chrombpnet_prep_nonpeaks:
#     input: "results/macs3_callpeak/{sample}-{genome}-.05/{sample}-{genome}-.05_peaks.narrowPeak"
#     output: "results/chrombpnet/{sample}-{genome}/chrombpnet_prep_nonpeaks.ok"
#     threads: 10
#     resources: 
#         slurm_partition='gpu', runtime=1440, mem_mb='96G', cpus_per_task=10, slurm_extra="--gres=gpu:1"
#     log: "logs/chrombpnet/prep_nonpeaks/{sample}-{genome}.log"
#     singularity: 'docker://kundajelab/chrombpnet:latest'
#     shell:
#         '''
#         (
#         chrombpnet prep nonpeaks -g resources/{wildcards.genome}.fa -p {input} -c resources/{wildcards.genome}.sizes.txt -fl resources/fold_0.json -br resources/{wildcards.genome}.blacklist.bed.gz -o results/chrombpnet/{wildcards.sample}-{wildcards.genome}
#         touch results/chrombpnet/{wildcards.sample}-{wildcards.genome}/chrombpnet_prep_nonpeaks.ok
#         ) &> {log}
#         '''


# rule chrombpnet_bias:
#     input: 
#         prep_nonpeaks_ok = "results/chrombpnet/{sample}-{genome}/chrombpnet_prep_nonpeaks.ok",
#         peaks= "results/macs3_callpeak/{sample}-{genome}-.05/{sample}-{genome}-.05_peaks.narrowPeak",
#         tagalign = "results/bedpe_tagalign/{sample}-{genome}.tagalign.gz"
        
#     output: "results/chrombpnet/{sample}-{genome}/chrombpnet_bias.ok"
#     threads: 10
#     resources: 
#         slurm_partition='gpu', runtime=1440, mem_mb='96G', cpus_per_task=10, slurm_extra="--gres=gpu:1"
#     log: "logs/chrombpnet/bias/{sample}-{genome}.log"
#     singularity: 'docker://kundajelab/chrombpnet:latest'
#     shell:
#         '''
#         (
#         rm -rf results/chrombpnet/{wildcards.sample}-{wildcards.genome}/bias_model
#         chrombpnet bias pipeline -d ATAC -itag {input.tagalign} -g resources/{wildcards.genome}.fa -c resources/{wildcards.genome}.sizes.txt -p {input.peaks} -n results/chrombpnet/{wildcards.sample}-{wildcards.genome}_negatives.bed -fl resources/fold_0.json -b .5 -o results/chrombpnet/{wildcards.sample}-{wildcards.genome}/bias_model
#         touch {output}
#         ) &> {log}
#         '''
    

    








        