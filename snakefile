configfile: "config.yaml"

rule all:
  input:
    expand("raw_qc/{sample}_{replicate}_fastqc.{extension}", sample=config["samples"], replicate=[1, 2], extension=["zip", "html"]),
    expand("Trimmed/{sample}_{replicate}{extension}.fastq.gz", sample=config["samples"], replicate=[1, 2], extension=["P", "U"]),
    "kallisto/Homo_sapiens.GRCh38.cdna.all.index",
    expand("kallisto/{condition}", condition=config["conditions"])


rule fastqc:
  input:
    rawread=expand("raw_data/{sample}_{replicate}.fastq.gz", sample=config["samples"], replicate=[1, 2])
  
  output:
    zip=expand("raw_qc/{sample}_{replicate}_fastqc.zip", sample=config["samples"], replicate=[1, 2]),
    internet=expand("raw_qc/{sample}_{replicate}_fastqc.html", sample=config["samples"], replicate=[1, 2])
  
  threads: 
    8
  
  params:
    path="raw_qc/"
  
  shell:
    "fastqc -t {threads} {input.rawread} -o {params.path}" 



rule trimmomatic:
  input:
    R1=expand("raw_data/{sample}_1.fastq.gz", sample=config["samples"]),
    R2=expand("raw_data/{sample}_2.fastq.gz", sample=config["samples"])
  
  output:
     RP=expand("Trimmed/{samples}_{replicate}P.fastq.gz", samples=config["samples"], replicate=[1, 2]),
     RU=expand("Trimmed/{samples}_{replicate}U.fastq.gz", samples=config["samples"], replicate=[1, 2])
  
  threads: 
    8 
  
  params:
    baseout=expand("Trimmed/{sample}.fastq.gz", sample=config["samples"]),
    log=expand("Trimmed/{sample}.log", sample=config["samples"]),
    adapter="adapter.fa"
  
  
  shell:
    "trimmomatic PE -threads {threads} -phred33 {input.R1} {input.R2} -baseout {params.baseout} \
    ILLUMINACLIP:{params.adapter}:2:30:10:2:keepBothReads LEADING:3 TRAILING:30 MINLEN:36 2>{params.log}"



rule kallisto_index:
  input:
    initial="kallisto/Homo_sapiens.GRCh38.cdna.all.fa"
    
  output:
    final="kallisto/Homo_sapiens.GRCh38.cdna.all.index"
    
  threads:
   8

  shell: "cp {input.initial} {output.final} | kallisto index -i {output.final} {input.initial}"     



rule kallisto_quant:
  input:
    read1=expand("Trimmed/{sample}_1P.fastq.gz", sample=config["samples"]),
    read2=expand("Trimmed/{sample}_2P.fastq.gz", sample=config["samples"]),
    index="kallisto/Homo_sapiens.GRCh38.cdna.all.index"
    
  output:
    expand("kallisto/{condition}", condition=config["conditions"])

  threads:
    8
  
  params:
    dir_out=expand("kallisto/{condition}", condition=config["conditions"]),

  shell: "kallisto quant -i {input.index} -o {params.dir_out} -t {threads} {input.read1} {input.read2}"
