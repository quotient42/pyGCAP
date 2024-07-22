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

1. `Python` >= 3.6
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

   - If you want to make a new conda environment for pygcap, follow the instructions below:

     ```
     conda create -n pygcap
     conda activate pygcap
     pip install pygcap
     conda install -c conda-forge ncbi-datasets-cli
     conda install -c conda-forge -c bioconda mmseqs2
     ```

---

## Usage

- pypi pygcap ([link](https://pypi.org/project/pygcap/))

  ```python
  # pip install pygcap
  pygcap [TAXON] [PROBE_FILE]
  ```

- input argument description

  ```python
  ### usage example
  pygcap Facklamia pygcap/data/probe_sample.tsv
  pygcap 66831 pygcap/data/probe_sample.tsv
  ```

  2.  `taxon` (both name and taxid are available)
  3.  path of `probe.tsv` ([sample file](https://github.com/jrim42/pyGCAP/blob/main/pygcap/data/probe_sample.tsv))

      - `Probe Name` (user defined)
      - `Prediction` (user defined)
      - `Accession` (UniProt entry)

### Options

1. `--working_dir` or `-w` (default: `.`): Specify the working directory path.

   ```
   pygcap [TAXON] [PROBE_FILE] —-working_dir or -w [PATH_OF_WORKING_DIRECTORY]
   ```

2. `--thread` or `-t` (default: `50`): Number of threads to use when running MMseqs2 and blastp. The number of threads can be adjusted automatically based on the CPU environment. It must be an integer greater than 0.

   ```
   pygcap [TAXON] [PROBE_FILE] —-thread or -t [NUMBER_OF_THREAD]
   ```

3. `--identity` of `-i` (default: `0.5`): The value of protein identity to be used in MMseqs2. It must be a value between 0 and 1.

   ```
   pygcap [TAXON] [PROBE_FILE] —-identity or -i [PROTEIN_IDENTITY]
   ```

4. `--skip` of `-s` (default: `none`): Specify steps to skip during the process. Multiple steps can be skipped by using this option multiple times.

   ```
   pygcap [TAXON] [PROBE_FILE] —-skip or -s [ARG]
   ```

   - `all`: Skip all the processes listed below.
   - `ncbi`: Skip downloading genome data from NCBI.
   - `mmseqs2`: Skip running MMseqs2.
   - `parsing`: Skip parsing genome data.
   - `uniprot`: Skip downloading probe data from UniProt.
   - `blastdb`: Skip running makeblastdb.

---

## (WIP)Output

- A directory with the following structure will be created in your `working directory` with the name of the `TAXON` provided as input.\

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
