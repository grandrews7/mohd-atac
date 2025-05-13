import pandas as pd
import os
import glob
import sys

data_dir = 'data/merged_data'
genomes = ['GRCh38']
qvals = ['.05']
samples = sorted([os.path.basename(x).replace("_R1.fastq.gz", "") for x in glob.glob(os.path.join(data_dir, '*_R1.fastq.gz'))])
print(samples)
samples = [x for x in samples if 'CKD_0003' in x]
samples = ['test']
print(f'There are {len(samples)} samples')
print(*samples)

reads = ["R1", "R2"]
print(f"There are {len(samples)} samples:\n\t", *samples)

rule all:
    input:
        #expand("resources/{genome}.fa{suffix}", genome = genomes, suffix = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"])
        expand("results/fastqc/raw/{sample}_{read}_fastqc.html", sample=samples, read=reads),
        expand("results/cutadapt/{sample}_R1.fastq.gz", sample=samples),
        # expand("results/fastqc/cutadapt/{sample}_{read}_fastqc.html", sample=samples, read=reads),
        expand("results/bowtie2_align_pe/{sample}-{genome}.bam", genome = genomes, sample=samples),
        expand("results/filter/{sample}-{genome}.bam", genome = genomes, sample=samples),
        expand("results/picard/{sample}-{genome}.bam", genome = genomes, sample=samples),
        expand("results/bedpe_tagalign/{sample}-{genome}.bedpe", genome = genomes, sample=samples),
        expand("results/macs3/{sample}-{genome}/{sample}-{genome}_peaks.narrowPeak", genome = genomes, sample=samples),
        expand("results/macs3_callpeak/{sample}-{genome}-{q}/{sample}-{genome}-{q}_peaks.narrowPeak", genome = genomes, sample=samples, q=qvals),
        expand("results/macs3/{sample}-{genome}/{sample}-{genome}.bigWig", genome = genomes, sample=samples),
        expand("results/qc/frag_len/{sample}-{genome}.png", genome = genomes, sample=samples), 
        expand("results/qc/frag_len/{sample}-{genome}.txt", genome = genomes, sample=samples),
        expand("results/qc/tss_enrichment/{sample}-{genome}.png", genome = genomes, sample=samples),
        expand("results/qc/tss_enrichment/{sample}-{genome}.txt", genome = genomes, sample=samples),
        expand("results/qc/frip_all/{sample}-{genome}-{q}.txt", genome = genomes, sample=samples, q=qvals),
        expand("results/chrombpnet/{sample}-{genome}/chrombpnet_prep_nonpeaks.ok", genome = genomes, sample=samples),
        expand("results/chrombpnet/{sample}-{genome}/chrombpnet_bias.ok", genome = genomes, sample=samples)
        # expand("results/pseudoreps/{sample}.pseudorep1.tagalign", sample=samples),
        # expand("results/pseudoreps/{sample}.pseudorep2.tagalign", sample=samples),
        # expand("results/macs3_pseudoreps/{sample}-{pseudorep}/{sample}-{pseudorep}_peaks.narrowPeak", sample=samples, pseudorep=["pseudorep1", "pseudorep2"]),
        # expand("results/idr/{sample}.IDR.bed", sample=samples),
        # expand("results/overlap_peaks/{sample}.overlap_peaks.bed", sample=samples),
        # expand("results/qc/tss_enrichment/{sample}.png", sample=samples),
        # expand("results/qc/tss_enrichment/{sample}.txt", sample=samples),
        # expand("results/qc/frag_len/{sample}.txt", sample=samples),
        # expand("results/qc/frag_len/{sample}.txt", sample=samples),
        # expand("results/qc/frip_all/{sample}.txt", sample=samples),
        # expand("results/qc/tss_enrichment/{sample}.png", sample=samples)
    
rule fastqc_raw:
    input: os.path.join(data_dir, "{sample}_{read}.fastq.gz")
    output: "results/fastqc/raw/{sample}_{read}_fastqc.html"
    conda:
        "MOHD-ATAC"
    threads: 8
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: "logs/fastqc_raw/{sample}-{read}.log"
    shell:
        """
        (
        fastqc -t {threads} -o results/fastqc/raw {input}
        ) &> {log}
        """

rule cutadapt:
    input: 
        R1 = os.path.join(data_dir, "{sample}_R1.fastq.gz"),
        R2 = os.path.join(data_dir, "{sample}_R2.fastq.gz")
    output:
        R1 = "results/cutadapt/{sample}_R1.fastq.gz",
        R2 = "results/cutadapt/{sample}_R2.fastq.gz"
    threads: 8
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    conda:
        "MOHD-ATAC"
    log:
        "logs/cutadapt/{sample}.log",
    
    shell:
        """
        (
        cutadapt -q 10 -m 15 -e 0.10 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -j {threads} -o {output.R1} -p {output.R2} {input.R1} {input.R2}
        ) &> {log}
        """
        
        
rule fastqc_cutadapt:
    input: "results/cutadapt/{sample}_{read}.fastq.gz"
    output: "results/fastqc/cutadapt/{sample}_{read}_fastqc.html"
    conda:
        "MOHD-ATAC"
    threads: 8
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: "logs/fastqc_trimmed/{sample}-{read}.log"
    shell:
        """
        (
        fastqc -t {threads} -o results/fastqc/cutadapt {input}
        ) &> {log}
        """

rule get_genome:
    output:
        fa = "resources/{genome}.fa",
        sizes = "resources/{genome}.sizes.txt",
    threads: 1
    log:
        "logs/get_genome/{genome}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        echo "Downloading genome fasta files"
        if [ {wildcards.genome} = 'GRCh38' ]
        then
            echo 'GRCh38'
            url='http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
        else
            echo 'T2T'
            url='http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/genome/t2t-chm13-v1.0.fa.gz'
        fi
        wget -q -O {output.fa}.gz $url
        gzip -d -c {output.fa}.gz > {output.fa}
        rm {output.fa}.gz
        samtools faidx {output.fa}
        cut -f1,2 {output.fa}.fai > {output.sizes}
        ) &> {log}
        """
        
rule bowtie2_build:
    input:
        fa = "resources/{genome}.fa",
    output:
        multiext("resources/{genome}.fa", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    threads: 24
    log:
        "logs/bowtie2_build/{genome}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        echo "Running bowtie2 build"

        bowtie2-build --threads {threads} {input} {input}
        ) &> {log}
        """
        
rule bowtie2_align_pe:
    input:
        R1 = "results/cutadapt/{sample}_R1.fastq.gz",
        R2 = "results/cutadapt/{sample}_R2.fastq.gz"

    output:
        "results/bowtie2_align_pe/{sample}-{genome}.bam"
    threads: 8 
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log:
        "logs/bowtie2_align_pe/{sample}-{genome}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        echo "Running bowtie2 align"
        bowtie2 -X 2000 --mm -k 1 \
        --threads {threads} \
        -x resources/{wildcards.genome}.fa \
        --rg-id {wildcards.sample}-{wildcards.genome} \
        --rg SM:{wildcards.sample}-{wildcards.genome} \
        -1 {input.R1} \
        -2 {input.R2} | \
        samtools view -1 -S /dev/stdin > /tmp/{wildcards.sample}-{wildcards.genome}.tmp.bam
        samtools sort -@ {threads} -T {wildcards.sample}-{wildcards.genome} -o {output} /tmp/{wildcards.sample}-{wildcards.genome}.tmp.bam
        samtools flagstat -@ {threads} {output} > {output}.flagstat
        rm /tmp/{wildcards.sample}-{wildcards.genome}.tmp.bam
        ) &> {log}
        """
        
rule filter:
    input: "results/bowtie2_align_pe/{sample}-{genome}.bam"
    output: "results/filter/{sample}-{genome}.bam"
    threads: 8 
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: 
        "logs/filter/{sample}-{genome}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        samtools view -F 1804 -f 2 -q 30 -u {input} | \
        samtools sort -n /dev/stdin -o /tmp/{wildcards.sample}-{wildcards.genome}.tmp.nmsrt.bam -T {wildcards.sample}-{wildcards.genome} -@ {threads}

        samtools view -h /tmp/{wildcards.sample}-{wildcards.genome}.tmp.nmsrt.bam |
        samtools fixmate -r /dev/stdin /tmp/{wildcards.sample}-{wildcards.genome}.tmp.fixmate.bam
        
        samtools view -F 1804 -f 2 -u /tmp/{wildcards.sample}-{wildcards.genome}.tmp.fixmate.bam | \
            samtools sort /dev/stdin -@ {threads} -o /tmp/{wildcards.sample}-{wildcards.genome}.tmp.coordsort.bam 

        samtools view -h /tmp/{wildcards.sample}-{wildcards.genome}.tmp.coordsort.bam | grep -v "chrM" | samtools view -b > {output}
        
        samtools index -@ {threads} {output} {output}.bai
        samtools flagstat -@ {threads} {output} > {output}.flagstat

        rm /tmp/{wildcards.sample}-{wildcards.genome}.tmp.nmsrt.bam
        rm /tmp/{wildcards.sample}-{wildcards.genome}.tmp.fixmate.bam
        rm /tmp/{wildcards.sample}-{wildcards.genome}.tmp.coordsort.bam
        ) &> {log}
        """

rule picard:
    input: "results/filter/{sample}-{genome}.bam"
    output: "results/picard/{sample}-{genome}.bam", "results/picard/{sample}-{genome}.txt"
    threads: 8 
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: 
        "logs/picard/{sample}-{genome}.log"
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
    input: "results/filter/{sample}-{genome}.bam"
    output: bedpe = "results/bedpe_tagalign/{sample}-{genome}.bedpe",
            tagalign = "results/bedpe_tagalign/{sample}-{genome}.tagalign.gz"
            
    threads: 8 
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: 
        "logs/bedpe_tagalign/{sample}-{genome}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        python scripts/bam_to_tagalign.py {input} {output.bedpe} results/bedpe_tagalign/{wildcards.sample}-{wildcards.genome}.tagalign --threads {threads}
        touch {output.tagalign}
        ) &> {log}
        """

rule macs3:
    input: "results/bedpe_tagalign/{sample}-{genome}.tagalign.gz"
    output:
        "results/macs3/{sample}-{genome}/{sample}-{genome}_peaks.narrowPeak", 
        "results/macs3/{sample}-{genome}/{sample}-{genome}_treat_pileup.bdg"
    threads: 8 
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    conda:
        "MOHD-ATAC"
    log: "logs/macs3/{sample}-{genome}.log"
    shell:
        """
        (
        macs3 callpeak -f BED -t {input} -n {wildcards.sample}-{wildcards.genome} --outdir results/macs3/{wildcards.sample}-{wildcards.genome} --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q .01 -B
        touch {output[0]}
        touch {output[1]}
        ) &> {log}
        """

rule macs3_callpeak:
    input: "results/bedpe_tagalign/{sample}-{genome}.tagalign.gz"
    output:
        "results/macs3_callpeak/{sample}-{genome}-{q}/{sample}-{genome}-{q}_peaks.narrowPeak"
    threads: 8 
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    conda:
        "MOHD-ATAC"
    log: "logs/macs3/{sample}-{genome}-{q}.log"
    shell:
        """
        (
        macs3 callpeak -f BED -t {input} -n {wildcards.sample}-{wildcards.genome}-{wildcards.q} --outdir results/macs3_callpeak/{wildcards.sample}-{wildcards.genome}-{wildcards.q} --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q {wildcards.q}
        touch {output}
        ) &> {log}
        """
        
        
rule bedGraphToBigWig:
    input: "results/macs3/{sample}-{genome}/{sample}-{genome}_treat_pileup.bdg"
    output: "results/macs3/{sample}-{genome}/{sample}-{genome}.bigWig"
    threads: 1
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='4G', cpus_per_task=1
    conda:
        "MOHD-ATAC"
    log: "logs/bedGraphToBigWig/{sample}-{genome}.log"
    shell:
        """
        (
        bedtools intersect -a {input} -b resources/{wildcards.genome}.sort.bed -sorted > results/macs3/{wildcards.sample}-{wildcards.genome}/{wildcards.sample}-{wildcards.genome}.bg
        bedGraphToBigWig results/macs3/{wildcards.sample}-{wildcards.genome}/{wildcards.sample}-{wildcards.genome}.bg resources/{wildcards.genome}.sizes.txt {output}
        ) &> {log}
        """

rule pseudoreps:
    input: "results/bedpe_tagalign/{sample}.bedpe"
    output: "results/pseudoreps/{sample}.pseudorep1.tagalign", "results/pseudoreps/{sample}.pseudorep2.tagalign"
    threads: 1
    conda:
        "MOHD-ATAC"
    log: "logs/pseudoreps/{sample}.log"
    shell:
        """
        (
        python scripts/split_bedpe_to_replicates.py {input} {output[0]} {output[1]}
        ) &> {log}
        """    

rule macs3_pseudoreps:
    input: "results/pseudoreps/{sample}.{pseudorep}.tagalign"
    output:
        "results/macs3_pseudoreps/{sample}-{pseudorep}/{sample}-{pseudorep}_peaks.narrowPeak"
    threads: 24
    conda:
        "MOHD-ATAC"
    log:"logs/macs3_pseudoreps/{sample}-{pseudorep}.log"
    shell:
        """
        (
        macs3 callpeak -f BED -t {input} -n {wildcards.sample}-{wildcards.pseudorep} --outdir results/macs3_pseudoreps/{wildcards.sample}-{wildcards.pseudorep} --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q .01 -B
        touch {output}
        ) &> {log}
        """


rule idr:
    input: "results/macs3_pseudoreps/{sample}-pseudorep1/{sample}-pseudorep1_peaks.narrowPeak", "results/macs3_pseudoreps/{sample}-pseudorep2/{sample}-pseudorep2_peaks.narrowPeak"
    output: "results/idr/{sample}.IDR.bed", "results/idr/{sample}.IDR.png"
    threads: 4
    conda:
        "idr"
    log: "logs/idr/{sample}.log"
    shell:
        """
        (
        sort -k8,8nr {input[0]} > {input[0]}.sorted
        sort -k8,8nr {input[1]} > {input[1]}.sorted
        idr --samples {input[0]}.sorted {input[1]}.sorted --input-file-type narrowPeak --rank p.value --output-file /tmp/{wildcards.sample}.IDR.unsorted.bed --plot
        cat /tmp/{wildcards.sample}.IDR.unsorted.bed | sort -k1,1 -k2,2n > {output[0]}
        mv /tmp/{wildcards.sample}.IDR.unsorted.bed.png {output[1]}
        ) &> {log}
        """

rule overlap_peaks:
    input: "results/macs3_pseudoreps/{sample}-pseudorep1/{sample}-pseudorep1_peaks.narrowPeak", 
            "results/macs3_pseudoreps/{sample}-pseudorep2/{sample}-pseudorep2_peaks.narrowPeak"
    output: "results/overlap_peaks/{sample}.overlap_peaks.bed"
    threads: 4
    conda:
        "MOHD-ATAC"
    log: "logs/overlap_peaks/{sample}.log"
    shell:
        """
        (
        bedtools intersect -a {input[0]} -b {input[1]} -u | sort -k1,1 -k2,2n > {output}
        ) &> {log}
        """

rule frag_len:
    input: "results/picard/{sample}-{genome}.bam"
    output: "results/qc/frag_len/{sample}-{genome}.png", "results/qc/frag_len/{sample}-{genome}.txt"
    threads: 8
    log: "logs/qc/frag_len/{sample}-{genome}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        python scripts/plot_fragment_length_distr.py {input} {threads} {output[0]} {output[1]}
        ) &> {log}
        """


rule tss_enrichment:
    input: "results/picard/{sample}-{genome}.bam"
    output: "results/qc/tss_enrichment/{sample}-{genome}.png", "results/qc/tss_enrichment/{sample}-{genome}.txt"
    threads: 8 
    resources:
            slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: "logs/qc/tss_enrichment/{sample}-{genome}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        python scripts/tss_enrichment.py {input} resources/ENCFF493CCB.bed.gz {output[0]} {output[1]}
        ) &> {log}
        """

rule frip_all:
    input: peaks = "results/macs3_callpeak/{sample}-{genome}-{q}/{sample}-{genome}-{q}_peaks.narrowPeak", reads = "results/bedpe_tagalign/{sample}-{genome}.tagalign.gz"
    output: "results/qc/frip_all/{sample}-{genome}-{q}.txt"
    threads: 8 
    resources:
        slurm_partition='12hours', runtime=240, mem_mb='32G', cpus_per_task=8
    log: "logs/qc/frip_all/{sample}-{genome}-{q}.log"
    conda:
        "MOHD-ATAC"
    shell:
        """
        (
        python scripts/calculate_frip.py {input.reads} {input.peaks} {output}
        ) &> {log}
        """ 

rule chrombpnet_prep_nonpeaks:
    input: "results/macs3_callpeak/{sample}-{genome}-.05/{sample}-{genome}-.05_peaks.narrowPeak"
    output: "results/chrombpnet/{sample}-{genome}/chrombpnet_prep_nonpeaks.ok"
    threads: 10
    resources: 
        slurm_partition='gpu', runtime=1440, mem_mb='96G', cpus_per_task=10, slurm_extra="--gres=gpu:1"
    log: "logs/chrombpnet/prep_nonpeaks/{sample}-{genome}.log"
    singularity: 'docker://kundajelab/chrombpnet:latest'
    shell:
        '''
        (
        chrombpnet prep nonpeaks -g resources/{wildcards.genome}.fa -p {input} -c resources/{wildcards.genome}.sizes.txt -fl resources/fold_0.json -br resources/{wildcards.genome}.blacklist.bed.gz -o results/chrombpnet/{wildcards.sample}-{wildcards.genome}
        touch results/chrombpnet/{wildcards.sample}-{wildcards.genome}/chrombpnet_prep_nonpeaks.ok
        ) &> {log}
        '''


rule chrombpnet_bias:
    input: 
        prep_nonpeaks_ok = "results/chrombpnet/{sample}-{genome}/chrombpnet_prep_nonpeaks.ok",
        peaks= "results/macs3_callpeak/{sample}-{genome}-.05/{sample}-{genome}-.05_peaks.narrowPeak",
        tagalign = "results/bedpe_tagalign/{sample}-{genome}.tagalign.gz"
        
    output: "results/chrombpnet/{sample}-{genome}/chrombpnet_bias.ok"
    threads: 10
    resources: 
        slurm_partition='gpu', runtime=1440, mem_mb='96G', cpus_per_task=10, slurm_extra="--gres=gpu:1"
    log: "logs/chrombpnet/bias/{sample}-{genome}.log"
    singularity: 'docker://kundajelab/chrombpnet:latest'
    shell:
        '''
        (
        rm -rf results/chrombpnet/{wildcards.sample}-{wildcards.genome}/bias_model
        chrombpnet bias pipeline -d ATAC -itag {input.tagalign} -g resources/{wildcards.genome}.fa -c resources/{wildcards.genome}.sizes.txt -p {input.peaks} -n results/chrombpnet/{wildcards.sample}-{wildcards.genome}_negatives.bed -fl resources/fold_0.json -b .5 -o results/chrombpnet/{wildcards.sample}-{wildcards.genome}/bias_model
        touch {output}
        ) &> {log}
        '''
    

    








        