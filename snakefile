configfile: "config.yaml"


ruleorder: fastqc > fastp > kallisto_index > kallisto_quant


rule all:
    input:
      expand("raw_qc/{sample}_{replicate}_fastqc.{extension}", sample=config["samples"], replicate=[1, 2], extension=["zip", "html"]),
      expand("trimmed/{sample}_{replicate}_trim.fastq.gz", sample=config["samples"], replicate=[1, 2]),
      "kallisto/Homo_sapiens.GRCh38.cdna.all.index",
      expand("kallisto/{condition}", condition=config["conditions"])


rule fastqc:
    input:
      rawread=expand("raw_data/{sample}_{replicate}.fastq.gz", sample=config["samples"], replicate=[1, 2])
    
    output:
      compress=expand("raw_qc/{sample}_{replicate}_fastqc.zip", sample=config["samples"], replicate=[1, 2]),
      net=expand("raw_qc/{sample}_{replicate}_fastqc.html", sample=config["samples"], replicate=[1, 2])
    
    threads: 
      8
    
    params:
      path="raw_qc/"
    
    shell:
      "fastqc -t {threads} {input.rawread} -o {params.path}" 


rule fastp:
    input:
      R1=expand("raw_data/{sample}_1.fastq.gz", sample=config["samples"]),
      R2=expand("raw_data/{sample}_2.fastq.gz", sample=config["samples"])
    
    output:
       RP=expand("trimmed/{samples}_1_trim.fastq.gz", samples=config["samples"]),
       RU=expand("trimmed/{samples}_2_trim.fastq.gz", samples=config["samples"])
     
    threads: 
      8 
    
    params:
      lenght=30
    
    shell:
      "fastp -i {input.R1} -I {input.R2} -o {output.RP} -O {output.RU} -l {params.lenght} --adapter_fasta adapter.fa"


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
      read1=expand("trimmed/{sample}_1_trim.fastq.gz", sample=config["samples"]),
      read2=expand("trimmed/{sample}_2_trim.fastq.gz", sample=config["samples"]),
      index="kallisto/Homo_sapiens.GRCh38.cdna.all.index"
      
    output:
      final=expand("kallisto/{condition}", condition=config["conditions"])

    threads:
      8

    shell: "kallisto quant -i {input.index} -o {output.final} -t {threads} {input.read1} {input.read2}"