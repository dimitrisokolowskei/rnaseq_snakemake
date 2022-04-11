configfile:"config.yaml"
print (config['samples'])

rule all:
    input:
      expand("raw_qc/{sample}_{replicate}_fastqc.{extension}", sample=config["samples"], replicate=[1, 2], extension=["zip", "html"]),
      expand("Trimmed/{sample}_{replicate}_trim.fastq.gz", sample=config["samples"], replicate=[1, 2]),
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


rule trimmomatic:
    input:
      R1=expand("raw_data/{sample}_1.fastq.gz", sample=config["samples"]),
      R2=expand("raw_data/{sample}_2.fastq.gz", sample=config["samples"])
    
    output:
       RP=expand("Trimmed/{samples}_1_trim.fastq.gz", samples=config["samples"]),
       RU=expand("Trimmed/{samples}_2_trim.fastq.gz", samples=config["samples"])
    
    threads: 
      8 
    
    params:
      lenght=36
    
    
    shell:
      "fastp -i {input.R1} -I {input.R2} -o {output.RP} -O {output.RU} -l {params.lenght} --detect_adapter_for_pe"


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
      read1=expand("Trimmed/{sample}_1_trim.fastq.gz", sample=config["samples"]),
      read2=expand("Trimmed/{sample}_2_trim.fastq.gz", sample=config["samples"]),
      index="kallisto/Homo_sapiens.GRCh38.cdna.all.index"
      
    output:
      expand("kallisto/{condition}", condition=config["conditions"])

    threads:
      8
    
    params:
      dir_out=expand("kallisto/{condition}", condition=config["conditions"]),

    shell: "kallisto quant -i {input.index} -o {params.dir_out} -t {threads} {input.read1} {input.read2}"
    
  
