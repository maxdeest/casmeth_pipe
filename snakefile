rule all:
    input:         
        "barcode10_methylation_frequencies_CYP19A1.tsv", 
        "barcode11_methylation_frequencies_CYP19A1.tsv", 
        "barcode12_methylation_frequencies_CYP19A1.tsv",
        
        
        #"barcode10_coverage_BDNF", "barcode05_coverage_TRPA1", "barcode10_coverage_CYP19A1",
        #"barcode06_coverage_BDNF", "barcode06_coverage_TRPA1", "barcode11_coverage_CYP19A1",
        #"barcode04_coverage_BDNF", "barcode04_coverage_TRPA1", "barcode12_coverage_CYP19A1",
        #"barcode07_coverage_BDNF", "barcode07_coverage_TRPA1", "barcode07_coverage_CYP19A1",
        
        
        

rule merge: 
    input:
        "/home/max/lab/ONT_seq/exp29/basecall/{sample}/"        #change
        
    output:
        "analysis/{sample}_analysis/{sample}.fastq"     
        
    shell:
        "cat {input}/*.fastq > {output}"

rule nano_index:
    input:
        "/home/max/lab/ONT_seq/exp29/fast5/",                   #change
        "analysis/{sample}_analysis/{sample}.fastq"
    output:
        "analysis/{sample}_analysis/{sample}.fastq.index", "analysis/{sample}_analysis/{sample}.fastq.index.fai",
        "analysis/{sample}_analysis/{sample}.fastq.index.gzi", "analysis/{sample}_analysis/{sample}.fastq.index.readdb"
    shell:
        "nanopolish index --directory {input}"

rule mini_align:
    input:
        ref="/home/max/Recommendref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        fa="analysis/{sample}_analysis/{sample}.fastq"
    output:
        "analysis/{sample}_analysis/{sample}.aln.sam"
    shell:
        "minimap2 -L -ax  map-ont {input.ref} {input.fa} > {output}" 

rule sam_to_bam:
    input:
        "analysis/{sample}_analysis/{sample}.aln.sam"
    output:
        "analysis/{sample}_analysis/{sample}.aln.bam"
    shell:
        "samtools view -S -b {input} > {output}"

rule sorting:
    input:
        "analysis/{sample}_analysis/{sample}.aln.bam"
    output:
        "analysis/{sample}_analysis/{sample}.aln.sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

rule sam_index:
    input:
        "analysis/{sample}_analysis/{sample}.aln.sorted.bam"
    output:
        "analysis/{sample}_analysis/{sample}.aln.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule call_meth:
    input:
        "analysis/{sample}_analysis/{sample}.aln.sorted.bam.bai",
        "analysis/{sample}_analysis/{sample}.fastq.index", "analysis/{sample}_analysis/{sample}.fastq.index.fai",
        "analysis/{sample}_analysis/{sample}.fastq.index.gzi", "analysis/{sample}_analysis/{sample}.fastq.index.readdb",
        r="analysis/{sample}_analysis/{sample}.fastq",
        b="analysis/{sample}_analysis/{sample}.aln.sorted.bam",
        g="/home/max/Recommendref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    params:
        w="chr11" #:27696900-27704000"
    
    output:
        "{sample}_methylation_calls_BDNF.tsv"
    shell:
        "nanopolish call-methylation -t 8 -r {input.r} -b {input.b} -g {input.g} -w {params.w} > {output}"

rule call_methtwo:
    input:
        "analysis/{sample}_analysis/{sample}.aln.sorted.bam.bai",
        "analysis/{sample}_analysis/{sample}.fastq.index", "analysis/{sample}_analysis/{sample}.fastq.index.fai",
        "analysis/{sample}_analysis/{sample}.fastq.index.gzi", "analysis/{sample}_analysis/{sample}.fastq.index.readdb",
        r="analysis/{sample}_analysis/{sample}.fastq",
        b="analysis/{sample}_analysis/{sample}.aln.sorted.bam",
        g="/home/max/Recommendref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    params:
        w="chr15" #:51277000-51292400"
    
    output:
        "{sample}_methylation_calls_CYP19A1.tsv"
    shell:
        "nanopolish call-methylation -t 8 -r {input.r} -b {input.b} -g {input.g} -w {params.w} > {output}"

rule call_meththree:
    input:
        "analysis/{sample}_analysis/{sample}.aln.sorted.bam.bai",
        "analysis/{sample}_analysis/{sample}.fastq.index", "analysis/{sample}_analysis/{sample}.fastq.index.fai",
        "analysis/{sample}_analysis/{sample}.fastq.index.gzi", "analysis/{sample}_analysis/{sample}.fastq.index.readdb",
        r="analysis/{sample}_analysis/{sample}.fastq",
        b="analysis/{sample}_analysis/{sample}.aln.sorted.bam",
        g="/home/max/Recommendref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    params:
        w="chr8"#:72043047-72089877"
    
    output:
        "{sample}_methylation_calls_TRPA1.tsv"
    shell:
        "nanopolish call-methylation -t 8 -r {input.r} -b {input.b} -g {input.g} -w {params.w} > {output}"
rule transform:
    input:
        '{sample}_methylation_calls_CYP19A1.tsv'
    output:
        '{sample}_methylation_frequencies_CYP19A1.tsv'
    shell:
        'python3 scripts/calculate_methylation_frequency.py {input} > {output}'

rule transformtwo:
    input:
        '{sample}_methylation_calls_BDNF.tsv'
    output:
        '{sample}_methylation_frequencies_BDNF.tsv'
    shell:
        'python3 scripts/calculate_methylation_frequency.py {input} > {output}'

rule transformthree:
    input:
        '{sample}_methylation_calls_TRPA1.tsv'
    output:
        '{sample}_methylation_frequencies_TRPA1.tsv'
    shell:
        'python3 scripts/calculate_methylation_frequency.py {input} > {output}'

rule coverageBDNF:
    input:
        "analysis/{sample}_analysis/{sample}.aln.sorted.bam"
    output:
        "{sample}_coverage_BDNF"
    shell:
        "samtools coverage -r chr11:27696900-27704000 -m {input} > {output}" 

rule coverageCYP19A1:
    input:
        "analysis/{sample}_analysis/{sample}.aln.sorted.bam"
    output:
        "{sample}_coverage_CYP19A1"
    shell:
        "samtools coverage -r chr15:51277000-51292400 -m {input} > {output}"

rule coverageTRPA1:
    input:
        "analysis/{sample}_analysis/{sample}.aln.sorted.bam"
    output:
        "{sample}_coverage_TRPA1"
    shell:
        "samtools coverage -r chr8:72043047-72089877 -m {input} > {output}"