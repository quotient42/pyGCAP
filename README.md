# lacto-dcw

- [introduction](#introduction)
- [program-flow](#program-flow)
- [pre-requirement](#pre-requirement)
- [usage](#usage)

---

### introduction

- Profiling _dcw_ genes from pan-genomes of Lactobacillales (LAB)

* this repository was built for the guiding 2023-2024 winter URAP student [@jrim42](https://github.com/jrim42)
* based on this paper - [Ancient origin and constrained evolution of the division and cell wall gene cluster in Bacteria](https://www.nature.com/articles/s41564-022-01257-y)
* Cell wall biosynthesis of LAB (eventually we can check and consider using related genes and proteins; for example, [Table 1](https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-13-S1-S9/tables/1) in this [reveiw paper](https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-13-S1-S9))
  <br>
  <img width=640 height = 480 src="https://media.springernature.com/lw685/springer-static/image/art%3A10.1186%2F1475-2859-13-S1-S9/MediaObjects/12934_2014_Article_1029_Fig2_HTML.jpg?as=webp">

---
### program-flow

<img width="801" alt="flowchart" src="https://github.com/logcossin/lacto-dcw/assets/90167645/db6e78d4-16b2-4a18-925d-c3dbb6fb3494">

---

### pre-requirement

  1. `Python`
  2. `conda` environment
      - `blast` [(bioconda blast package)](https://anaconda.org/bioconda/blast)
        ```
        conda install bioconda::blast
        conda install bioconda/label/cf201901::blast
        ```
      - `datasets` & `dataformat` from NCBI [(conda-forge ncbi-datasets-cli package)](https://anaconda.org/conda-forge/ncbi-datasets-cli)
        
          ```
          conda install conda-forge::ncbi-datasets-cli
          ```
      - `MMseqs2` ([MMseqs2 github](https://github.com/soedinglab/MMseqs2))
        ```
        conda install -c conda-forge -c bioconda mmseqs2
        ```

### usage
  ```
  pgc-finder <working_dir> <TAXON> <target_path>
  ```
  
  - input data description
     1. working directory
     2. taxon (both name and taxid are available)
     3. path of `target.tsv`
         - target name (user defined)
         - protein name (user defined)
         - protein entry (UniProt) 

