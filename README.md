# PGCfinder
PGCfinder: A Python Package for Probe-based Gene Cluster Finding in Large Microbial Genome Database

- [introduction](#introduction)
- [pipeline-flow](#pipeline-flow)
- [pre-requirement](#pre-requirement)
- [usage](#usage)
- [example](#example)

---

### introduction

Bacterial gene clusters provide insights into metabolism and evolution, and facilitate biotechnological applications. We developed PGCfinder, a Python package for probe-based gene cluster discovery. This pipeline uses sequence search and analysis tools and public databases (e.g. BLAST,  MMSeqs2,  UniProt, and NCBI) to predict potential gene clusters by user-provided probe genes. We tested the pipeline with the division and cell wall (dcw) gene cluster, crucial for cell division and peptidoglycan biosynthesis. 

To evaluate PGCfinder, we used 17 major dcw genes defined by Megrian et al. [1] as a probe set to search for gene clusters in 696 Lactobacillales genomes. The results were integrated to provide detailed information on gene content,  gene order, and types of clusters. While PGCfinder examined the completeness of the gene clusters, it could also suggest novel taxa-specific accessory genes related to dcw clusters in Lactobacillales genomes. The package will be freely available on the Python Package Index, Bioconda, and GitHub.

[1] Megrian, D., et al. [Ancient origin and constrained evolution of the division and cell wall gene cluster in Bacteria](https://www.nature.com/articles/s41564-022-01257-y). Nat Microbiol 7, 2114–2127 (2022).
  
---
### pipeline-flow

<p align="center">
  <img width="801" alt="flowchart" src="https://github.com/logcossin/lacto-dcw/assets/90167645/db6e78d4-16b2-4a18-925d-c3dbb6fb3494">
</p>

---

### pre-requirement

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

### usage
  ```
  pip install pgc-finder
  pgc-finder <working_dir> <TAXON> <target_path>
  ```
  - pypi pgc-finder [(link)](https://pypi.org/project/pgc-finder/)
  - input data description
     1. working directory
     2. taxon (both name and taxid are available)
     3. path of `target.tsv`
         - target name (user defined)
         - protein name (user defined)
         - protein entry (UniProt) 

---

### example

- Profiling _dcw_ genes from pan-genomes of Lactobacillales (LAB)
* Cell wall biosynthesis of LAB (eventually we can check and consider using related genes and proteins; for example, [Table 1](https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-13-S1-S9/tables/1) in this [reveiw paper](https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-13-S1-S9))
  <br>
  <img width=640 height = 480 src="https://media.springernature.com/lw685/springer-static/image/art%3A10.1186%2F1475-2859-13-S1-S9/MediaObjects/12934_2014_Article_1029_Fig2_HTML.jpg?as=webp">

