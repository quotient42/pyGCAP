## (work in process)

# lacto-dcw

- [introduction](#introduction)
- [pre-requirement](#pre-requirement)
- [tutorial](#tutorial)
- [result description](#result-description)
- [example](#example)

---

### introduction

- Profiling _dcw_ genes from pan-genomes of Lactobacillales (LAB)

* this repository was built for the guiding 2023-2024 winter URAP student [@jrim42](https://github.com/jrim42)
* based on this paper - [Ancient origin and constrained evolution of the division and cell wall gene cluster in Bacteria](https://www.nature.com/articles/s41564-022-01257-y)
* Cell wall biosynthesis of LAB (eventually we can check and consider using related genes and proteins; for example, [Table 1](https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-13-S1-S9/tables/1) in this [reveiw paper](https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-13-S1-S9))
  <br>
  <img width=640 height = 480 src="https://media.springernature.com/lw685/springer-static/image/art%3A10.1186%2F1475-2859-13-S1-S9/MediaObjects/12934_2014_Article_1029_Fig2_HTML.jpg?as=webp">

---

### pre-requirement

- 필요한 프로그램

  1. `Python`
  2. `conda` environment
  3. `blast` (CLI)
  4. `datasets` & `dataformat` from NCBI
  5. `MMseqs2` ([github](https://github.com/soedinglab/MMseqs2))

- 필요한 data

  1. target protein data
     - `target.tsv`
     - `target.fasta`
  2. ncbi dataset -> target organism의 `genomic.gff`, `genomic.gbff`, `protein.faa`

     ```
     # Get tools used for data download
     wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
     wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat

     # Change the file mode
     chmod a+x ./datasets # file should appear green when viewed with 'ls' command
     chmod a+x ./dataformat
     ```

     ```
     ./datasets download genome taxon "lactobacillales" --reference --include protein,gff3 --filename {OUTPUT_FILENAME}.zip

     unzip {OUTPUT_FILENAME}.zip
     ```

  3. ncbi dataformat -> `assembly_report.tsv`
     ```
     ./dataformat tsv genome --inputfile {YOUR_PATH}/ncbi_dataset/data/assembly_data_report.jsonl --fields accession,organism-name,assmstats-total-sequence-len,checkm-completeness,checkm-contamination,annotinfo-featcount-gene-total > assembly_report.tsv
     ```
  4. MMseqs2 -> `result_clutster.tsv`

     ```
     cat {YOUR_PATH}/ncbi_dataset/data/*/protein.faa > {OUTPUT_FILENAME}.faa

     {YOUR_PATH}/mmseqs/bin/mmseqs easy-cluster {OUTPUT_FILENAME}.faa result tmp --min-seq-id 0.5 --threads 55
     ```

---

### tutorial

1. 아래 명령어 실행하여 project dir 지정하고 필요한 정보를 넣기
   ```
   python3 src/main/init.py
   ```
   ```
   project root dir
   ├── data          assembly_report.tsv, target.fasta, result_clutster.tsv
   ├── input         all GCF_* directories
   ├── output
   │   ├── tsv
   │   ├── img
   │   └── genus
   └── seqlib
   ```
   - project의 초기 설정은 `src/main`의 `project_info`에서 확인가능
2. 아래 명령어 실행하여 target contents/order 구하기
   ```
   python3 src/main/main.py
   ```

---

### result description

1. `output/tsv`
1. `output/img`
1. `output/genus`

---

### exmaple
