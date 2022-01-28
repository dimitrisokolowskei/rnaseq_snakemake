# RNA-Seq Analysis Workflow Automation Using Snakemake
## Introduction
RNA sequencing (RNA-Seq) employs high-troughtput sequencing technologies to unravel the transcriptomic profile of living organisms. The raw sequencing *output* *dataset* requires a series of elegant and sofisticated pieces of bioinformatics *softwares* to gather usefull and insightful information on biological data. Here I'll describe a basic RNA-Seq *pipeline* using snakemake workflow language for automation.
## *Software* Requirements 
The *softwares* required to run this analysis are contained in a conda environment. If you don't have anaconda installed in your machine, a basic tutorial can be foud [here](https://www.digitalocean.com/community/tutorials/how-to-install-the-anaconda-python-distribution-on-ubuntu-20-04). Otherwise, you can simply import an enviroment that I've already made using:
```sh
conda env create --file enviroment.yaml
conda activate rnaseq
```
This command will use the ```enviroment.yaml``` file and create a ```rnaseq``` enviroment containing all the programs you'll need. 