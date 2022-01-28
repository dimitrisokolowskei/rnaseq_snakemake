# RNA-Seq Analysis Workflow Automation Using Snakemake
## Introduction
RNA sequencing (RNA-Seq) employs high-troughtput sequencing technologies to unravel the transcriptomic profile of living organisms. The raw sequencing *dataset* requires a series of elegant and sofisticated pieces of bioinformatics *softwares* to gather usefull and insightful information on biological data. Here I'll describe a basic RNA-Seq *pipeline* using snakemake workflow language for *script* automation.
## *Software* Requirements 
The *softwares* required to run this analysis are contained in a conda environment. If you don't have anaconda installed in your machine, a basic tutorial can be foud [here](https://www.digitalocean.com/community/tutorials/how-to-install-the-anaconda-python-distribution-on-ubuntu-20-04). Otherwise, you can simply import an enviroment that I've already made using:
```sh
conda env create --file environment.yaml # Create rnaseq environment
conda activate rnaseq # enters the recently created environment
```
This command will use the ```enviroment.yaml``` file and create a enviroment called ```rnaseq``` containing all the programs you'll need. 

## Data Download 
The RNA-Seq data that we're going to use is described in [Sousa et al., 2019](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5967-8). In this paper, the autors describes DEG of T cells stimulated by different versions of OKT3, an anti-CD3 antibody. This data is pubicly available in SRA with the study code ```SRP139131```. For matters of simplification we'll only use the control and the T cells treated with OKT3 duplicates data. However, feel free to use the entire dataset. First let's create our working sub(directories) using: 
```sh
mkdir -p rnaseq/{Trimmed,metadata,kallisto,raw_data,raw_qc}
```
Now, inside ```raw_data``` download our work data from SRA using the following command:
```sh
cat SRR_Acc_List.txt | parallel "fastq-dump --gzip --split-files {}"
```
The reference genome also need to be collected for latter mapping process. Inside ```kallisto``` directory use:
```sh
curl -O http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz | gunzip
```



