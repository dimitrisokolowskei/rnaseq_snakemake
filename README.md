# RNA-Seq Analysis Workflow Automation Using Snakemake
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.github.io)
## Introduction
RNA sequencing (RNA-Seq) employs high-troughtput sequencing technologies to unravel the transcriptomic profile of living organisms. The raw sequencing *dataset* requires a series of elegant and sofisticated pieces of bioinformatics *softwares* to gather usefull and insightful information on biological data. Here I'll describe a basic RNA-Seq *pipeline* using snakemake workflow language for automation.
## *Software* Requirements 
The *softwares* required to run this analysis are contained in a conda environment. If you don't have anaconda installed in your machine, a basic tutorial can be foud [here](https://www.digitalocean.com/community/tutorials/how-to-install-the-anaconda-python-distribution-on-ubuntu-20-04). Otherwise, you can simply import an enviroment that I've already made using:
```sh
conda env create --file environment.yaml # Create rnaseq environment
conda activate rnaseq # enters the recently created environment
```
This command will use the ```enviroment.yaml``` file and create a enviroment called ```rnaseq``` containing all the programs you'll need. 

## Data Download 
The RNA-Seq data that we're going to use is described in [Sousa et al., 2019](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5967-8). In this paper, the autors describes DEG of T cells stimulated by different versions of OKT3, an anti-CD3 antibody. This data is pubicly available in SRA with the study code ```SRP139131```. For matters of simplification we'll only use the control and the T cells treated with OKT3 data. However, feel free to use the entire dataset or any other piece of data that may be of your interest.

Inside ```rnaseq``` directory, use the following command to create our (sub)directories: 
```sh
mkdir -p {Trimmed,metadata,kallisto,raw_data,raw_qc}
```
Now, inside ```raw_data``` directory, download our work data from SRA using the following command:
```sh
cat SRR_Acc_List.txt | parallel "fastq-dump --gzip --split-files {}"
```
The reference genome also need to be collected for latter mapping process. Inside ```kallisto``` directory use:
```sh
curl -O http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz && gunzip *.gz
```
## Snakemake Run
To run our analysis we'll need to execute the ```snakefile``` containing all the RNA-Seq steps. In your terminal, use:
```sh
snakemake --core 8 all
```
```snakemake``` command will execute our ```snakefile``` file. The number of cores after ```--core``` will depend on your machine or server CPU capability. 
## Warnings 
- Before executing ```snakefile``` ensure that all the data and directories needed are as described earlier.
- If you try to use your own personal data, be aware that the ```config.yaml``` file needs to be modified. You should only insert the files name/id, instead of something like ```{sample_name}_1.fastq.gz```. 
- In the ```config.yaml``` file the number of ```samples``` **needs** to be equal to the number of ```conditions```. And every ```samples``` **must** have a ```conditions``` with a diferent name. This happens duo to kallisto pseudoaligner, where *output* files will be sent to a specific directory (```conditions```). 
- If you are interested to only run a specific step of the ```snakefile```, just specify the name of the **rule** you are interest to run, instead of ```all```. Example:
 ```sh
snakemake --core 8 fastqc
```
