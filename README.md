# pyGCAP: a (py)thon (G)ene (C)luster (A)nnotation & (P)rofiling

A Python Package for Probe-based Gene Cluster Finding in Large Microbial Genome Database

- [Introduction](#Introduction)
- [Pipeline-flow](#Pipeline-flow)
- [Pre-requirement](#Pre-requirement)
- [Usage](#Usage)

---

## Introduction

Bacterial gene clusters provide insights into metabolism and evolution, and facilitate biotechnological applications. We developed pyGCAP, a Python package for probe-based gene cluster discovery. This pipeline uses sequence search and analysis tools and public databases (e.g. BLAST, MMSeqs2, UniProt, and NCBI) to predict potential gene clusters by user-provided probe genes. We tested the pipeline with the division and cell wall (dcw) gene cluster, crucial for cell division and peptidoglycan biosynthesis.

To evaluate pyGCAP, we used 17 major dcw genes defined by Megrian et al. [1] as a probe set to search for gene clusters in 716 Lactobacillales genomes. The results were integrated to provide detailed information on gene content, gene order, and types of clusters. While PGCfinder examined the completeness of the gene clusters, it could also suggest novel taxa-specific accessory genes related to dcw clusters in Lactobacillales genomes. The package will be freely available on the Python Package Index, Bioconda, and GitHub.

[1] Megrian, D., et al. [Ancient origin and constrained evolution of the division and cell wall gene cluster in Bacteria](https://www.nature.com/articles/s41564-022-01257-y). Nat Microbiol 7, 2114–2127 (2022).

---

## Pipeline-flow

<p align="center">
  <img width="1000" alt="flowchart" src="https://github.com/jrim42/pyGCAP/assets/90167645/a39af39e-7961-4e21-b2ab-e1a3c86b1f4a">
</p>

---

## Pre-requirement

1. `Python`
2. `conda` environment

   - `blast` ([bioconda blast package](https://anaconda.org/bioconda/blast))
     ```
     conda install bioconda::blast
     conda install bioconda/label/cf201901::blast
     ```
   - `datasets` & `dataformat` from NCBI ([conda-forge ncbi-datasets-cli package](https://anaconda.org/conda-forge/ncbi-datasets-cli))

     ```
     conda install conda-forge::ncbi-datasets-cli
     ```

   - `MMseqs2` ([MMseqs2 github](https://github.com/soedinglab/MMseqs2))
     ```
     conda install -c conda-forge -c bioconda mmseqs2
     ```

---

## Usage

- pypi pygcap ([link](https://pypi.org/project/pygcap/))

  ```python
  pip install pygcap
  conda activate ncbi_datasets
  pygcap [WORKING_DIRECTORY] [TAXON] [PROBE_FILE]
  ```

- input argument description

  ```python
  ### usage example
  pygcap . Facklamia pygcap/data/probe_sample.tsv
  pygcap . 66831 pygcap/data/probe_sample.tsv
  ```

  1.  `working directory`
  2.  `taxon` (both name and taxid are available)
  3.  path of `probe.tsv` ([sample file](https://github.com/jrim42/pyGCAP/blob/main/pygcap/data/probe_sample.tsv))

      - `Probe Name` (user defined)
      - `Prediction` (user defined)
      - `Accession` (UniProt entry)

### Options

- `--skip`: Specify steps to skip during the process. Multiple steps can be skipped by using this option multiple times.

  - `all`: to skip all processes below
  - `ncbi`: to skip downloading genome data from `NCBI`
  - `mmseqs2`: to skip running `MMseqs2`
  - `parsing`: to skip parsing genome data
  - `uniprot`: to skip downloading probe data from `Uniprot`
  - `blastdb`: to skip running `makeblastdb`

  ```
  pygcap [WORKING_DIRECTORY] [TAXON] [PROBE_FILE] —-skip or -s [ARG]
  ```

---

## (WIP)Output

- A directory with the following structure will be created in your `working directory` with the name of the `TAXON` provided as input.
  ```
  📦 [TAXON_NAME]
  ├─ data
  │  ├─ assembly_report.tsv
  │  ├─ metadata_target.tsv
  │  └─ ...
  ├─ input
  │  ├─ [GENUS_01]
  │  ├─ [GENUS_02]
  │  └─ ...
  ├─ output
  │  ├─ genus
  │  ├─ img
  │  └─ tsv
  └─ seqlib
     ├─ blast_output.tsv
     ├─ seqlib.tsv
     └─ ...
  ```

---

## (WIP) example

- Profiling _dcw_ genes from pan-genomes of Lactobacillales (LAB)
