configfile: "config.yaml"


WORK_DATA = config["WORK_DATA"]
WORK_QC = config["WORK_QC"]
WORK_TRIM = config["WORK_TRIM"]
WORK_KALL = config["WORK_KALL"]


rule all:
    input:
      expand(WORK_QC + "{sample}_{replicate}_fastqc.{extension}", sample=config["samples"], replicate=[1, 2], extension=["zip", "html"]),
      expand(WORK_TRIM + "{sample}_{replicate}_trim.fastq.gz", sample=config["samples"], replicate=[1, 2]),
      WORK_KALL + "Homo_sapiens.GRCh38.cdna.all.index",
      expand(WORK_KALL + "quant_result_{condition}", condition=config["conditions"])
      


      
rule fastqc:
    input:
      WORK_DATA + "{sample}_{replicate}.fastq.gz"
    
    output:
      WORK_QC + "{sample}_{replicate}_fastqc.{extension}",
      
    priority: 50

    threads: 
      8
    
    params:
      path = WORK_QC
    
    shell:
      "fastqc -t {threads} {input} -o {params.path}" 



rule fastp:
    input:
      R1 = WORK_DATA + "{sample}_1.fastq.gz",
      R2 = WORK_DATA + "{sample}_2.fastq.gz"
    
    output:
       RP = WORK_TRIM + "{sample}_1_trim.fastq.gz",
       RU = WORK_TRIM + "{sample}_2_trim.fastq.gz"
       
    priority: 40
     
    threads: 
      8 
    
    params:
      lenght=30
    
    shell:
      "fastp -i {input.R1} -I {input.R2} -o {output.RP} -O {output.RU} -l {params.lenght} --adapter_fasta adapter.fa"


rule kallisto_index:
    input:
      WORK_KALL + "Homo_sapiens.GRCh38.cdna.all.fa"
      
    output:   
      WORK_KALL + "Homo_sapiens.GRCh38.cdna.all.index"
      
    threads:
     8

    shell: 
      "cp {input} {output} | kallisto index -i {output} {input}"      



rule kallisto_quant:
    input:
      fq1 = WORK_TRIM + "{sample}_1_trim.fastq.gz",
      fq2 = WORK_TRIM + "{sample}_2_trim.fastq.gz",
      idx = WORK_KALL + "Homo_sapiens.GRCh38.cdna.all.fa.index"
    
    output:
      WORK_KALL + "quant_result_{condition}"

    threads: 
       8
    
    shell:
      "kallisto quant -i {input.idx} -t {threads} -o {output} {input.fq1} {input.fq2}"