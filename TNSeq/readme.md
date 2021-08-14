# DNA-seq workflow

This workflow is mainly based on the commercial software Sentieon. The workflow performs bioinformatics pipeline for Tumor-Normal analysis recommended in the Broad institute Somatic short variant discovery (SNVs + Indels). It uses Docker containers making installation trivial and results highly reproducible. All docker images used in the workflow are stored on XDP, and the workflow can be run on XDP directly. The following figure illustrates such a typical bioinformatics pipeline.
![img.png](img.png)

## Table of Contents

- [Summary](#Summary)
- [Install](#install)
- [Tools](#tools)
- [Usage](#usage)
- [To do list](#todolist)
- [Maintainers](#maintainers)
- [License](#license)

## Summary

The workflow has the following main steps:
* bwa_mem: alignment to reference genome
* get_metrics and plot_metrics: Calculate data metrics and plot the metrics
* DeDup: remove/mark duplicate reads
* realign: perform in-del realignment
* recalibration: Base quality score recalibration (BQSR).
* Somatic variant discovery, with the following stages:
    1. (Optional) Estimate the cross-sample contamination and tumor segmentation.
    2. (Optional) Estimate any possible orientation bias present in the sequencing.
    3. Somatic candidate variant calling on the two individual BAM files: this step identifies the potential sites where the cancer genome data displays somatic variations relative to the normal genome, and calculates genotypes at that site.
    4. Filter the variants.
* Annotate variants with the SnpEff software

## Install

This workflow is developed to be run on XDP without any installation.

## Tools
* sentieon
* snpEff:5.0ef

## Usage

Log in XDP using Google Chrome browser. Go to the section "Mine/App" and create the workflow by providing two files: the argument specification file and WDL workflow file. After creating the workflow application, workflow instance can be launched in project section.

### Files for creating workflow App on XDP
- Argument specification file: *args.detail.xlsx*
- WDL workflow file: pipeline described in [WDL](https://github.com/openwdl/wdl)

### Argument detail
- please refer to **args.detail.xlsx**

### All input files of the workflow
In this bioinformatics pipeline you will need the following inputs:
- The FASTA file containing the nucleotide sequence of the reference genome corresponding to the sample you will analyze.
- Two sets of FASTQ files containing the nucleotide sequence of the sample to be analyzed, one for the tumor sample and one for the matched normal sample. These files contain the raw reads from the DNA sequencing. The software supports inputting FASTQ files compressed using GZIP. The software only supports files containing quality scores in Sanger format (Phred+33).
- You can also include in the pipeline the following optional inputs that will help the algorithms detect artifacts and remove false positives:
    * Panel of normal VCF: list of common errors that appear as variants from multiple unrelated normal samples. The contents of this file will be used to identify variants that are more likely to be germline variants, and filter them as such.
    * Population resource VCF: list of population allele specific frequencies that will be used for filtering possible germline variants and to annotate the results.
- snpEff database, you may download corresponding annotation database according to the manual of snpEff 



## To do list
- support multiple paired samples
- support tumor only mode
- support CNV detection
- support VEP annotation

## Maintainers

Email: *danny.gu@basebit.ai*

## License

[MIT © Richard McRichface.](../LICENSE)