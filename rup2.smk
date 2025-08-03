import os
from itertools import product


(PLANTS, )=glob_wildcards('{plants}/reference/genome.fa')
(PLANTS2, SAMPLESPE, PS, )=glob_wildcards('{plants}/reads/{samples}_{s}.fq.gz')

def filter_combinator(combinator):
    def filtered_combinator(*args, **kwargs):
        print("filter files...")
        for wc_comb in combinator(*args, **kwargs):
            filetocheck_1 = "%s/reads/%s_1.fq.gz" % (wc_comb[0][1], wc_comb[1][1])
            filetocheck_2 = "%s/reads/%s_2.fq.gz" % (wc_comb[0][1], wc_comb[1][1])
            filetocheck_s = "%s/reads/%s_s.fq.gz" % (wc_comb[0][1], wc_comb[1][1])
            
            if(os.path.isfile(filetocheck_1) and os.path.isfile(filetocheck_2)):
                yield wc_comb
            if(os.path.isfile(filetocheck_s)):
                yield wc_comb

    return filtered_combinator

filtered_product = filter_combinator(product)



rule all:
    input:
        expand('{plant}/{plant}.report.pdf', plant=PLANTS)



rule index:
    input:
        expand('{plant}/results/index/salmon/index.complete', plant=PLANTS),
        expand('{plant}/results/index/salmon_rrna/index.complete', plant=PLANTS),
        expand('{plant}/results/index/salmon_quantiles/index.complete', plant=PLANTS)



rule trim:
    input:
        expand('{plant}/results/trimmed/{sample}.json',
               filtered_product, plant=PLANTS, sample=SAMPLESPE)
           


rule quant:
    input:
        expand('{plant}/results/salmon/{sample}/quant.sf',
               filtered_product, plant=PLANTS, sample=SAMPLESPE)
    output:
        '{plant}/results/salmon/.done'
    shell:
        "touch {output}"



rule quant_small:
    input:
        expand('{plant}/esults/salmon_{ref}/{sample}/quant.sf',
               filtered_product, plant=PLANTS, sample=SAMPLESPE, ref=['rrna', 'quantiles'])
    output:
        '{plant}/results/salmon_{ref}/.done'
    shell:
        "touch {output}"



rule index_full:
    input:
        t="{plant}/results/index/salmon/gentrome.fa",
        decoys="{plant}/results/index/salmon/decoys.txt"
    output:
        "{plant}/results/index/salmon/index.complete"
    threads:
        24
    resources:
        mem_mb=256000,
        runtime=360,
        nodes=1,
        cpus_per_task=24
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check index -p {threads} \
           -i {wildcards.plant}/results/index/salmon \
                  -t {input.t}  --decoys {input.decoys} --keepDuplicates && \
           touch {output}
           """



rule index_small:
    input:
        t="{plant}/reference/{ref}.fa",
    output:
        "{plant}/results/index/salmon_{ref}/index.complete"
    threads:
        8
    resources:
        mem_mb=16000,
        runtime=360,
        nodes=1,
        cpus_per_task=8
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check index -p {threads} \
                  -i {wildcards.plant}/results/index/salmon_{wildcards.ref} \
                  -t {input.t} --keepDuplicates && \
           touch {output}
           """



rule gentrome:
    input:
        genome="{plant}/reference/genome.fa",
        fai="{plant}/reference/genome.fa.fai",
        transcripts="{plant}/reference/transcripts.fa"
    output:
        gentrome="{plant}/results/index/salmon/gentrome.fa",
        decoys="{plant}/results/index/salmon/decoys.txt"
    shell: """
           cat {input.transcripts} {input.genome} > {output.gentrome}
           cut -f 1 {input.fai} > {output.decoys}
           """



rule transcripts:
    input:
        genome="{plant}/reference/genome.fa",
        annotation="{plant}/reference/annotation.gtf"
    output:
        "{plant}/reference/transcripts.fa"
    conda: "envs/gffread.yaml"
    shell: """
           gffread {input.annotation} -g {input.genome} -w {output} 
           """



rule quantiles:
    input:
        fasta="{plant}/reference/transcripts.fa"
    output:
        quantiles="{plant}/reference/quantiles.fa"
    script:
        "scripts/split10.py"



rule rrna:
    input:
        genome="{plant}/reference/genome.fa",
        annotation="{plant}/reference/rRNA.gtf"
    output:
        "{plant}/reference/rrna.fa"
    conda: "envs/gffread.yaml"
    shell: """
           gffread {input.annotation} -g {input.genome} -w {output}
           """



rule find_rrnas:
    input: "{plant}/reference/genome.fa"
    output: "{plant}/reference/barrnap.out"
    conda: "envs/barrnap.yaml"
    threads: 8
    shell: """barrnap --kingdom euk --threads {threads} {input} > {output}"""


           
rule rrna_gtf:
    input: "{plant}/reference/barrnap.out"
    output: "{plant}/reference/rRNA.gtf"
    script: "scripts/barrnap2gtf.py"


            
rule trim_pe:
    input:
        r1="{plant}/reads/{sample}_1.fq.gz",
        r2="{plant}/reads/{sample}_2.fq.gz"
    output:
        r1pe="{plant}/results/trimmed/{sample}_R1.fastq.gz",
        r2pe="{plant}/results/trimmed/{sample}_R2.fastq.gz",
        r1se="{plant}/results/trimmed/{sample}_s1.fastq.gz",
        r2se="{plant}/results/trimmed/{sample}_s2.fastq.gz",
        json="{plant}/results/trimmed/{sample}.json"
    threads:
        8
    conda: "envs/fastp.yaml"
    shell: """
           fastp -i {input.r1} -I {input.r2} \
                 -o {output.r1pe} -O {output.r2pe} \
                 --unpaired1 {output.r1se} --unpaired2 {output.r2se} \
                 -j {output.json} -l 50 -h /dev/null
           """



rule trim_se:
    input:
        rs="{plant}/reads/{sample}_s.fq.gz"
    output:
        rse="{plant}/results/trimmed/{sample}_S.fastq.gz",
        json="{plant}/results/trimmed/{sample}.json"
    threads:
        8
    conda: "envs/fastp.yaml"
    shell: """
           fastp -i {input.rs}  \
                 -o {output.rse}  \
                 -j {output.json} -l 50 -h /dev/null
           """



rule salmon_pe:
    input:
        r1pe="{plant}/results/trimmed/{sample}_R1.fastq.gz",
        r2pe="{plant}/results/trimmed/{sample}_R2.fastq.gz",
        index="{plant}/results/index/salmon/index.complete"
    output:
        '{plant}/results/salmon/{sample}/quant.sf'
    threads:
        8
    resources:
        mem_mb=64000,
        runtime=360,
        nodes=1,
        cpus_per_task=8
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check quant -l A --numGibbsSamples 30 \
		  --gcBias --validateMappings \
		  --allowDovetail \
	   	  -i {wildcards.plant}/results/index/salmon \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon/{wildcards.sample}\
		  -1 {input.r1pe} -2 {input.r2pe}
           """



rule salmon_se:
    input:
        rse="{plant}/results/trimmed/{sample}_S.fastq.gz",
        index="{plant}/results/index/salmon/index.complete"
    output:
        '{plant}/results/salmon/{sample}/quant.sf'
    threads:
        8
    resources:
        mem_mb=64000,
        runtime=360,
        nodes=1,
        cpus_per_task=8
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check quant -l A --numGibbsSamples 30 \
		  --gcBias --validateMappings \
		  --allowDovetail \
	   	  -i {wildcards.plant}/results/index/salmon \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon/{wildcards.sample}\
		  -r {input.rse}
           """



rule salmon_pe_rrna:
    input:
        r1pe="{plant}/results/trimmed/{sample}_R1.fastq.gz",
        r2pe="{plant}/results/trimmed/{sample}_R2.fastq.gz",
        index="{plant}/results/index/salmon_rrna/index.complete"
    output:
        '{plant}/results/salmon_rrna/{sample}/quant.sf'
    threads:
        8
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check quant -l A \
	   	  -i {wildcards.plant}/results/index/salmon_rrna \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon_rrna/{wildcards.sample}\
		  -1 {input.r1pe} -2 {input.r2pe}
           """




rule salmon_se_rrna:
    input:
        rse="{plant}/results/trimmed/{sample}_S.fastq.gz",
        index="{plant}/results/index/salmon_rrna/index.complete"
    output:
        '{plant}/results/salmon_rrna/{sample}/quant.sf'
    threads:
        8
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check quant -l A \
	   	  -i {wildcards.plant}/results/index/salmon_rrna \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon_rrna/{wildcards.sample}\
		  -r {input.rse} 
           """



rule salmon_pe_quantiles:
    input:
        r1pe="{plant}/results/trimmed/{sample}_R1.fastq.gz",
        r2pe="{plant}/results/trimmed/{sample}_R2.fastq.gz",
        index="{plant}/results/index/salmon_quantiles/index.complete"
    output:
        '{plant}/results/salmon_quantiles/{sample}/quant.sf'
    threads:
        8
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check quant -l A \
	   	  -i {wildcards.plant}/results/index/salmon_quantiles \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon_quantiles/{wildcards.sample}\
		  -r {input.r1pe} {input.r2pe}
           """



rule salmon_se_quantiles:
    input:
        rse="{plant}/results/trimmed/{sample}_S.fastq.gz",
        index="{plant}/results/index/salmon_quantiles/index.complete"
    output:
        '{plant}/results/salmon_quantiles/{sample}/quant.sf'
    threads:
        8
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check quant -l A \
	   	  -i {wildcards.plant}/results/index/salmon_quantiles \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon_quantiles/{wildcards.sample}\
		  -r {input.rse}
           """



rule fai:
    input: "{genome}"
    output: "{genome}.fai"
    conda: "envs/samtools.yaml"
    shell: "samtools faidx {input}"



rule report:
    input:
        '{plant}/reference/samples.tsv',
        expand('{plant}/results/trimmed/{sample}.json',
               filtered_product, plant=PLANTS, sample=SAMPLESPE),
        expand('{plant}/results/salmon/{sample}/quant.sf',
               filtered_product, plant=PLANTS, sample=SAMPLESPE),
        expand('{plant}/results/salmon_{ref}/{sample}/quant.sf',
               filtered_product, plant=PLANTS, sample=SAMPLESPE, ref=['rrna', 'quantiles'])
    output:
        report='{plant}/{plant}.report.pdf'
    conda: "envs/R.yaml"
    params: reads=100000
    script: "scripts/rup2.R"



rule init_conda:
    input:
        ".envs_samtools",
        ".envs_salmon",
        ".envs_fastp",
        ".envs_gffread",
        ".envs_R"
    output:
        ".conda"
    shell: """touch {output}"""



rule init_salmon:
    output: ".envs_salmon"
    conda: "envs/salmon.yaml"
    shell: """salmon --version > {output}"""



rule init_fastp:
    output: ".envs_fastp"
    conda: "envs/fastp.yaml"
    shell: """fastp --version > {output}"""



rule init_gffread:
    output: ".envs_gffread"
    conda: "envs/gffread.yaml"
    shell: """gffread 2> {output} || echo"""



rule init_samtools:
    output: ".envs_samtools"
    conda: "envs/samtools.yaml"
    shell: """samtools --version > {output}"""



rule init_R:
    output: ".envs_R"
    conda: "envs/R.yaml"
    shell: """R --version > {output}"""