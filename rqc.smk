import os
from itertools import product
from pathlib import Path

localrules: report, gentrome, gtf, transcripts, quantiles, rrna, rrna_gtf, fai, plant_bam, init_conda, init_salmon, init_fastp, init_gffread, init_samtools, init_R


(PLANTS, GROUPS, )=glob_wildcards('{plants}/reference/samples{groups}.tsv')
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
                case "gtf":
                    yield("{plant}/results/scallop/{sample}/scallop.gtf").format(plant=p[0], sample=p[1])
                case "qc_untrimmed":
                    if(os.path.isfile(filetocheck_s)):
                        yield("{plant}/results/falco/untrimmed/{sample}_s.fq.gz").format(plant=p[0], sample=p[1])
                    else:
                        yield("{plant}/results/falco/untrimmed/{sample}_1.fq.gz").format(plant=p[0], sample=p[1])
                        yield("{plant}/results/falco/untrimmed/{sample}_2.fq.gz").format(plant=p[0], sample=p[1])
                case "qc_trimmed":
                    if(os.path.isfile(filetocheck_s)):
                        yield("{plant}/results/falco/trimmed/{sample}_S.fastq.gz").format(plant=p[0], sample=p[1])
                    else:
                        yield("{plant}/results/falco/trimmed/{sample}_R1.fastq.gz").format(plant=p[0], sample=p[1])
                        yield("{plant}/results/falco/trimmed/{sample}_R2.fastq.gz").format(plant=p[0], sample=p[1])
                case _:
                    yield("{plant}/results/{type}/{sample}/quant.sf").format(plant=p[0], type=type, sample=p[1])


                    
def get_sample_report_files():
    for p in zip(PLANTS, GROUPS):
        samples_file = "%s/reference/samples%s.tsv" % (p[0], p[1])

        if(os.path.isfile(samples_file)):
           yield("{plant}/report/samples{group}.report.html").format(plant=p[0], group=p[1])

    for p in PLANTS2:
        samples_file = "%s/reference/samples.tsv" % (p)
        if(os.path.isfile(samples_file)):
            yield("{plant}/report/samples.report.html").format(plant=p)



MAX_MEM_MB = config["max_mem_mb"]

def use_large(wc):
    genome = Path(f"{wc.plant}/reference/genome.fa")

    required_mb = ( genome.stat().st_size / 1024**3 * 13 * 1024 )

    return required_mb <= MAX_MEM_MB



rule report:
    input: lambda wildcards: get_sample_report_files()



rule report_subproject:
    input:
        samples='{plant}/reference/samples{group}.tsv',
        trimmed=lambda wildcards: get_inputs(wildcards.plant, "trimmed"),
        salmon=lambda wildcards: get_inputs(wildcards.plant, "salmon"),
        salmon_rrna=lambda wildcards: get_inputs(wildcards.plant, "salmon_rrna"),
        salmon_quantiles=lambda wildcards: get_inputs(wildcards.plant, "salmon_quantiles"),
        qc_untrimmed=lambda wildcards: get_inputs(wildcards.plant, "qc_untrimmed"),
        qc_trimmed=lambda wildcards: get_inputs(wildcards.plant, "qc_trimmed")
    output:
        report='{plant}/report/samples{group}.report.html'
    wildcard_constraints:
        group = ".*"
    conda: "envs/R.yaml"
    threads: 1
    benchmark: "{plant}/benchmark/samples{group}.report_html.txt"
    script: "scripts/rqc.R"

           

rule index_full:
    input:
        t="{plant}/results/index/salmon/gentrome.fa",
        decoys="{plant}/results/index/salmon/decoys.txt"
    output:
        "{plant}/results/index/salmon/index.complete"
    threads:
        24
    resources:
        mem_mb=MAX_MEM_MB,
        runtime=360,
        nodes=1,
        cpus_per_task=24
    conda: "envs/salmon.yaml"
    benchmark: "{plant}/benchmark/index_full.txt"
    shell: """
           DECOY="--decoys {input.decoys}"
           if [ ! -s {input.decoys} ] ; then
              DECOY=""
           fi
           salmon --no-version-check index -p {threads} \
                  -i {wildcards.plant}/results/index/salmon \
                  -t {input.t} $DECOY --keepDuplicates && \
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
        mem_mb=MAX_MEM_MB,
        runtime=360,
        nodes=1,
        cpus_per_task=8
    conda: "envs/salmon.yaml"
    benchmark: "{plant}/benchmark/index_{ref}.txt"
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
    params:
        mode=lambda wc: "large" if use_large(wc) else "small"
    benchmark: "{plant}/benchmark/gentrome.txt"
    shell: """
           if [ "{params.mode}" = "large" ] ; then
                cat {input.transcripts} {input.genome} > {output.gentrome}
                cut -f 1 {input.fai} > {output.decoys}
           else
                cat {input.transcripts} > {output.gentrome}
                touch {output.decoys}
                DECOY=$( basename {input.transcripts} )/decoys.fa
                if [ -e $DECOY ] ; then
                     cat $DECOY >> {output.gentrome}
                     grep "^>" $DECOY | sed "s/>//; s/ .*//" > {output.decoys}
                fi
           fi
           """



rule gtf:
    input: "{plant}/reference/annotation.gff3"
    output: "{plant}/reference/annotation.gtf"
    conda: "envs/gffread.yaml"
    shell: """
           gffread {input} -T -o {output} 
           """



rule transcripts:
    input:
        genome="{plant}/reference/genome.fa",
        annotation="{plant}/reference/transcripts.gtf"
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
    benchmark: "{plant}/benchmark/barrnap.txt"
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
    params: min_length = config.get("mininal_read_length", 30)
    benchmark: "{plant}/benchmark/trim_pe.{sample}.txt"
    shell: """
           fastp -i {input.r1} -I {input.r2} \
                 -o {output.r1pe} -O {output.r2pe} \
                 --unpaired1 {output.r1se} --unpaired2 {output.r2se} \
                 --thread {threads} \
                 -j {output.json} -l {params.min_length} -h /dev/null
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
    params: min_length = config.get("mininal_read_length", 30)
    benchmark: "{plant}/benchmark/trim_se.{sample}.txt"
    shell: """
           fastp -i {input.rs}  \
                 -o {output.rse}  \
                 --thread {threads} \
                 -j {output.json} -l {params.min_length} -h /dev/null
           """



def get_salmon_pe_input(wc):
    yield "{plant}/results/trimmed/{sample}_R1.fastq.gz"
    yield "{plant}/results/trimmed/{sample}_R2.fastq.gz"
    yield "{plant}/results/index/salmon/index.complete"

    if config["use_alignment"] == "true":
        yield "{plant}/results/bam/{sample}/STARAligned.sorted.bam"
    else:
        yield "{plant}/results/index/salmon/index.complete"

def get_salmon_se_input(wc):
    yield "{plant}/results/trimmed/{sample}_S.fastq.gz"
    yield "{plant}/results/index/salmon/index.complete"

    if config["use_alignment"] == "true":
        yield "{plant}/results/bam/{sample}/STARAligned.sorted.bam"
    else:
        yield "{plant}/results/index/salmon/index.complete"


rule salmon_pe:
    input: get_salmon_pe_input
    output:
        '{plant}/results/salmon/{sample}/quant.sf'
    threads:
        16 
    resources:
        mem_mb=MAX_MEM_MB,
        runtime=360,
        nodes=1,
        cpus_per_task=16
    conda: "envs/salmon.yaml"
    benchmark: "{plant}/benchmark/salmon_pe.{sample}.txt"
    params: use_alignment=config["use_alignment"]
    shell: """
           if [ "{params.use_alignment}" == "true" ] ; then
                INPUT="-a {input[3]}"
           else
                INPUT="-1 {input[0]} -2 {input[1]} -i {wildcards.plant}/results/index/salmon"
           fi
    
           salmon --no-version-check quant -l A --numGibbsSamples 30 \
		  --gcBias --validateMappings --minAssignedFrags 0 \
		  --allowDovetail \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon/{wildcards.sample}\
                  $INPUT
           """



rule salmon_se:
    input: get_salmon_se_input
    output:
        '{plant}/results/salmon/{sample}/quant.sf'
    threads:
        16
    resources:
        mem_mb=MAX_MEM_MB,
        runtime=360,
        nodes=1,
        cpus_per_task=16
    conda: "envs/salmon.yaml"
    benchmark: "{plant}/benchmark/salmon_se.{sample}.txt"
    params: use_alignment=config["use_alignment"]
    shell: """
           if [ "{params.use_alignment}" == "true" ] ; then
                INPUT="-a {input[2]}"
           else
                INPUT="-r {input[0]} -i {wildcards.plant}/results/index/salmon"
           fi

           salmon --no-version-check quant -l A --numGibbsSamples 30 \
		  --gcBias --validateMappings --minAssignedFrags 0 \
		  --allowDovetail \
	   	  -i {wildcards.plant}/results/index/salmon \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon/{wildcards.sample}\
		  $INPUT
           """



rule salmon_small_base:
    threads:
        16
    conda: "envs/salmon.yaml"
    shell: """
           salmon --no-version-check quant -l A --minAssignedFrags 0 \
	   	  -i {wildcards.plant}/results/index/salmon_{wildcards.type} \
		  -p {threads} \
		  -o {wildcards.plant}/results/salmon_{wildcards.type}/{wildcards.sample} \
		  -r {input.reads}
           """


use rule salmon_small_base as salmon_small_se with:
    input:
        reads="{plant}/results/trimmed/{sample}_S.fastq.gz",
        index="{plant}/results/index/salmon_{type}/index.complete"
    output:
        '{plant}/results/salmon_{type}/{sample}/quant.sf'
    benchmark: "{plant}/benchmark/salmon_se_{type}.{sample}.txt"


use rule salmon_small_base as salmon_small_pe with:
    input:
        reads=["{plant}/results/trimmed/{sample}_R1.fastq.gz",
               "{plant}/results/trimmed/{sample}_R2.fastq.gz"],
        index="{plant}/results/index/salmon_{type}/index.complete"
    output:
        '{plant}/results/salmon_{type}/{sample}/quant.sf'
    benchmark: "{plant}/benchmark/salmon_se_{type}.{sample}.txt"



rule fai:
    input: "{genome}"
    output: "{genome}.fai"
    conda: "envs/samtools.yaml"
    shell: "samtools faidx {input}"



rule falco_base:
    threads: 1
    conda: "envs/falco.yaml"
    shell: """
           falco {input} -skip-report -q \
                 -D {output.data} -S {output.summary}
           """


use rule falco_base as falco_untrimmed with:
    input: "{plant}/reads/{name}"
    output:
        data="{plant}/results/falco/untrimmed/{name}",
        summary="{plant}/results/falco/untrimmed/{name}.summary"
    benchmark: "{plant}/benchmark/falco_untrimmed.{name}.txt"

    
use rule falco_base as falco_trimmed with:
    input: "{plant}/results/trimmed/{name}"
    output:
        data="{plant}/results/falco/trimmed/{name}",
        summary="{plant}/results/falco/trimmed/{name}.summary"
    benchmark: "{plant}/benchmark/falco_trimmed.{name}.txt"



if config["use_assembly"] == "true":
    rule copy_scallop_gtf:
        input: "{plant}/results/scallop.gtf"
        output: "{plant}/reference/transcripts.gtf"
        shell: """ cp {input} {output} """
else:
    rule copy_user_gtf:
        input: "{plant}/reference/annotation.gtf"
        output: "{plant}/reference/transcripts.gtf"
        shell: """ cp {input} {output} """


               

#### CHECK CONDA ######

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




#### MAPPING / GUIDED ASSEMBLY #####

rule all_plants_bam:
    input: expand("{plant}/.bam_done", plant=PLANTS2)
    output: ".bam_done"
    shell: """ touch {output} """


rule plant_bam:
    input: lambda wildcards: get_inputs(wildcards.plant, "bam")
    output: "{plant}/.bam_done"
    shell: """ touch {output} """

  

rule assemblies:
    input: expand("{plant}/results/scallop.gtf", plant=PLANTS2)
    output: ".assembly_done"
    shell: """ touch {output} """






rule assembly:
    input:
        gtf=lambda wildcards: get_inputs(wildcards.plant, "gtf"),
    output: "{plant}/results/scallop.gtf"
    threads: 12
    benchmark: "{plant}/benchmark/stringtie_merge.txt"
    resources:
        mem_mb=15000,
        runtime=360,
        nodes=1,
        cpus_per_task=12
    conda: "envs/scallop.yaml"
    shell: """
           USEREF=""
           [ -e "{wildcards.plant}/reference/annotation.gtf" ] && USEREF="-G {wildcards.plant}/reference/annotation.gtf"
           stringtie --merge -o {output} $USEREF -p 12 --fr --conservative {input.gtf}
           """



rule star_map_base:
    threads: 24
    benchmark: "{plant}/benchmark/star_map.{sample}.txt"
    resources:
        mem_mb=96000,
        runtime=360,
        nodes=1,
        cpus_per_task=24
    conda: "envs/STAR.yaml"
    shell: """
                #### --twopassMode Basic 
           [ -e {wildcards.plant}/results/bam/{wildcards.sample}_tmp ] && \
               rm -fr {wildcards.plant}/results/bam/{wildcards.sample}_tmp

           STAR --genomeDir {input.index} \
                --outTmpDir {wildcards.plant}/results/bam/{wildcards.sample}_tmp \
                --readFilesIn {input.r1pe} {input.r2pe} \
                --outSAMunmapped Within \
                --outFilterType BySJout \
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

            rm -fr {wildcards.plant}/results/bam/{wildcards.sample}_tmp
            """

use rule star_map_base as star_map_pe with:
    input: 
        r1pe="{plant}/results/trimmed/{sample}_R1.fastq.gz",
        r2pe="{plant}/results/trimmed/{sample}_R2.fastq.gz",
        index="{plant}/results/index/STAR"
    output:
        "{plant}/results/bam/{sample}/STARAligned.out.bam"

use rule star_map_base as star_map_se with:
    input: 
        r1pe="{plant}/results/trimmed/{sample}_S.fastq.gz",
        r2pe=[],
        index="{plant}/results/index/STAR"
    output:
        "{plant}/results/bam/{sample}/STARAligned.out.bam"


rule sort_bam:
    input: "{plant}/results/bam/{sample}/STARAligned.out.bam"
    output: "{plant}/results/bam/{sample}/STARAligned.sorted.bam"
    threads: 8
    resources:
        mem_mb=64000,
        runtime=360,
        nodes=1,
        cpus_per_task=8
    benchmark: "{plant}/benchmark/sort_bam.{sample}.txt"
    conda: "envs/STAR.yaml"
    shell: """
            samtools sort -m 5G -@ {threads} -o {output} {input}

            samtools index {output}
            """



rule scallop:
    input: "{plant}/results/bam/{sample}/STARAligned.sorted.bam"
    output: "{plant}/results/scallop/{sample}/scallop.gtf"
    threads: 1
    benchmark: "{plant}/benchmark/scallop.{sample}.txt"
    resources:
        mem_mb=32000,
        runtime=360,
        nodes=1,
        cpus_per_task=1
    conda: "envs/scallop.yaml"
    shell: """
            scallop2 -i {input} -o {output} --threads {threads} --verbose 0
           """
           


rule star_index:
    input:
        genome="{plant}/reference/genome.fa",
    output: directory("{plant}/results/index/STAR")
    threads: 24
    benchmark: "{plant}/benchmark/star_index.txt"
    resources:
        mem_mb=256000,
        runtime=360,
        nodes=1,
        cpus_per_task=24
    conda: "envs/STAR.yaml"
    shell: """
      [ -e {wildcards.plant}/STARtmp ] && rm -fr {wildcards.plant}/STARtmp

      USEREF=""
      [ -e "{wildcards.plant}/reference/annotation.gtf" ] && USEREF="--sjdbGTFfile {wildcards.plant}/reference/annotation.gtf"
    
	   STAR --genomeSAindexNbases 11 \
                --runThreadN {threads} \
                --runMode genomeGenerate \
                --genomeFastaFiles {input.genome} \
                --sjdbGTFfile $USEREF \
                --sjdbOverhang 257 \
                --limitGenomeGenerateRAM 256000000000 \
                --outTmpDir {wildcards.plant}/STARtmp \
                --genomeDir {output}

      rm -fr {wildcards.plant}/STARtmp
           """
