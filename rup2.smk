import os
from itertools import product


(PLANTS, )=glob_wildcards('{plants}/reference/genome.fa')
(PLANTS2, SAMPLES, PE, )=glob_wildcards('{plants}/reads/{samples}_{s}.fq.gz')



def get_inputs(plant, type):
    for p in zip(PLANTS2, SAMPLES):
            filetocheck_1 = "%s/reads/%s_1.fq.gz" % (p[0], p[1])
            filetocheck_2 = "%s/reads/%s_2.fq.gz" % (p[0], p[1])
            filetocheck_s = "%s/reads/%s_s.fq.gz" % (p[0], p[1])

            if(p[0] == plant and ((os.path.isfile(filetocheck_1) and os.path.isfile(filetocheck_2)) or os.path.isfile(filetocheck_s))):
                match type:
                    case "trimmed":
                        yield("{plant}/results/trimmed/{sample}.json").format(plant=p[0], sample=p[1])
                    case "bam":
                        yield("{plant}/results/bam/{sample}/STARAligned.sorted.bam").format(plant=p[0], sample=p[1])
                    case _:
                        yield("{plant}/results/{type}/{sample}/quant.sf").format(plant=p[0], type=type, sample=p[1])



rule report:
    input:
        expand('{plant}/{plant}.report.pdf', plant=PLANTS)



rule index:
    input:
        expand('{plant}/results/index/salmon/index.complete', plant=PLANTS),
        expand('{plant}/results/index/salmon_rrna/index.complete', plant=PLANTS),
        expand('{plant}/results/index/salmon_quantiles/index.complete', plant=PLANTS)



rule trim:
    input:
        expand('{plant}/results/trimmed/{sample}.json', zip, plant=PLANTS2, sample=SAMPLES)
               


rule quant:
    input:
        expand('{plant}/results/salmon/{sample}/quant.sf', zip, plant=PLANTS2, sample=SAMPLES)



rule quant_rrna:
    input:
        expand('{plant}/results/salmon_rrna/{sample}/quant.sf', zip, plant=PLANTS2, sample=SAMPLES)



rule quant_quantiles:
    input:
        expand('{plant}/results/salmon_quantiles/{sample}/quant.sf', zip, plant=PLANTS2, sample=SAMPLES)



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


rule star_map:
    input: 
        r1pe="{plant}/results/trimmed/{sample}_R1.fastq.gz",
        r2pe="{plant}/results/trimmed/{sample}_R2.fastq.gz",
        index="{plant}/results/index/STAR"
    output:
        "{plant}/results/bam/{sample}/STARAligned.out.bam"
    threads: 24
    resources:
        mem_mb=96000,
        runtime=360,
        nodes=1,
        cpus_per_task=24
    conda: "envs/STAR.yaml"
    shell: """
           STAR --genomeDir {input.index} \
                --readFilesIn {input.r1pe} {input.r2pe} \
                --outSAMunmapped Within \
                --outFilterType BySJout \
                --twopassMode Basic \
                --outSAMattributes NH HI AS NM MD \
                --outFilterMultimapNmax 20 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverReadLmax 0.04 \
                --alignIntronMin 20 \
                --alignIntronMax 5000000 \
                --alignMatesGapMax 1000000 \
                --readFilesCommand zcat \
                --runThreadN {threads} \
                --outSAMtype BAM Unsorted \
                --outSAMstrandField intronMotif \
                --outSAMheaderHD @HD VN:1.4 SO:coordinate \
                --outFileNamePrefix {wildcards.plant}/results/bam/{wildcards.sample}/STAR
            """



rule sort_bam:
    input: "{plant}/results/bam/{sample}/STARAligned.out.bam"
    output: "{plant}/results/bam/{sample}/STARAligned.sorted.bam"
    threads: 8
    resources:
        mem_mb=64000,
        runtime=360,
        nodes=1,
        cpus_per_task=8
    conda: "envs/STAR.yaml"
    shell: """
            samtools sort -m 5G -@ {threads} -o {output} {input}

            samtools index {output}
            """



rule star_index:
    input:
        genome="{plant}/reference/genome.fa",
        gff="{plant}/reference/annotation.gtf"
    output: directory("{plant}/results/index/STAR")
    threads: 24
    resources:
        mem_mb=256000,
        runtime=360,
        nodes=1,
        cpus_per_task=24
    conda: "envs/STAR.yaml"
    shell: """
	   STAR --genomeSAindexNbases 11 \
                --runThreadN {threads} \
                --runMode genomeGenerate \
                --genomeFastaFiles {input.genome} \
                --sjdbGTFfile {input.gff} \
                --sjdbOverhang 257 \
                --limitGenomeGenerateRAM 256000000000 \
                --genomeDir {output}
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
                  -t {input.t} && \
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
        16 
    resources:
        mem_mb=64000,
        runtime=360,
        nodes=1,
        cpus_per_task=16
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
        16
    resources:
        mem_mb=64000,
        runtime=360,
        nodes=1,
        cpus_per_task=16
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
        16
    resources:
        mem_mb=32000,
        runtime=360,
        nodes=1,
        cpus_per_task=16
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check quant -l A \
	   	  -i {wildcards.plant}/results/index/salmon_rrna \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon_rrna/{wildcards.sample} \
		  -1 {input.r1pe} -2 {input.r2pe}
           """




rule salmon_se_rrna:
    input:
        rse="{plant}/results/trimmed/{sample}_S.fastq.gz",
        index="{plant}/results/index/salmon_rrna/index.complete"
    output:
        '{plant}/results/salmon_rrna/{sample}/quant.sf'
    threads:
        16
    resources:
        mem_mb=32000,
        runtime=360,
        nodes=1,
        cpus_per_task=16
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
        16
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check quant -l A \
	   	  -i {wildcards.plant}/results/index/salmon_quantiles \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon_quantiles/{wildcards.sample} \
		  -r {input.r1pe} {input.r2pe}
           """



rule salmon_se_quantiles:
    input:
        rse="{plant}/results/trimmed/{sample}_S.fastq.gz",
        index="{plant}/results/index/salmon_quantiles/index.complete"
    output:
        '{plant}/results/salmon_quantiles/{sample}/quant.sf'
    threads:
        16
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check quant -l A \
	   	  -i {wildcards.plant}/results/index/salmon_quantiles \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon_quantiles/{wildcards.sample} \
		  -r {input.rse}
           """



rule fai:
    input: "{genome}"
    output: "{genome}.fai"
    conda: "envs/samtools.yaml"
    shell: "samtools faidx {input}"



rule plant_bam:
    input: lambda wildcards: get_inputs(wildcards.plant, "bam")
    output: "{plant}/.bam_done"
    shell: """ touch {output} """



rule report_pdf:
    input:
        '{plant}/reference/samples.tsv',
        lambda wildcards: get_inputs(wildcards.plant, "trimmed"),
        lambda wildcards: get_inputs(wildcards.plant, "salmon"),
        lambda wildcards: get_inputs(wildcards.plant, "salmon_rrna"),
        lambda wildcards: get_inputs(wildcards.plant, "salmon_quantiles")
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
