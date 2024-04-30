(work in process)

# Profiling _dcw_ genes from pan-genomes of Lactobacillales

- this repository was built for the guiding 2023-2024 winter URAP student [@jrim42](https://github.com/jrim42)
- based on this paper - [Ancient origin and constrained evolution of the division and cell wall gene cluster in Bacteria](https://www.nature.com/articles/s41564-022-01257-y)
- Cell wall biosynthesis of LAB (eventually we can check and consider using related genes and proteins; for example, [Table 1](https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-13-S1-S9/tables/1) in this [reveiw paper](https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-13-S1-S9))
  <br>
  <img width=640 height = 480 src="https://media.springernature.com/lw685/springer-static/image/art%3A10.1186%2F1475-2859-13-S1-S9/MediaObjects/12934_2014_Article_1029_Fig2_HTML.jpg?as=webp">

---

### Pre-requirement

1. Download target genomes (Refseq Lactobacillales genomes)

```
# Get tools used for data download
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat

# Change the file mode
chmod a+x ./datasets # file should appear green when viewed with 'ls' command
chmod a+x ./dataformat
```

2. probe genes (_dcw_ genes and MurJ)
3. search tools
4. data analysis tools
5. visualization tools

### Working environment

1. `panpyro` server
2. `conda` environment
3. `Python` script and `R` package

### Algorithmic procedure & checking points

1. Download reference genomes of lactic acid bacteria (Lactobacillales, LAB) from NCBI RefSeq database with NCBI Datasets [command-line tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
   - 697 LAB genomes available in RefSeq as of 2024.01.04
   - Downloaded data are found in this folder `/home/jsrim/data`
   - **!!!!CHECK!!!** the quality of genomes (sizes/completeness/genes of genomes; number of species/taxa; 697 Refeq represents ###,### genomes of same taxa..) - make a summary table of target data <---- what Jisu should do

```
# Download LAB genome data
# Change the bold name to your prefered name, the default filename is ncbi_dataset.zip

./datasets download genome taxon "lactobacillales" --reference --include protein,gff3 --filename OUTPUT_FILENAME.zip

# Unzip downloaded data (should create folder named 'ncbi_dataset')

unzip OUTPUT_FILENAME.zip

# Get metadata of downloaded genomes

./dataformat tsv genome --inputfile ncbi_dataset/data/assembly_data_report.jsonl --fields accession,organism-name,assmstats-total-sequence-len,checkm-completeness,checkm-contamination,annotinfo-featcount-gene-total > assembly_report.tsv
```

2. Cluster LAB protein sequences with [MMseqs2](https://github.com/soedinglab/MMseqs2) with the referenece _dcw_ genes from [UniProt](https://www.uniprot.org/) or supplementary table
   - Cluster results are found in this folder `panpyro/data/...` [link]()
   - **CHECK** how probe genes were selected and parameters for clustering/model building

```
# Combine protein sequences to one file
# Change the bold name to your prefered name

cat /home/jsrim/data/ncbi_dataset/data/*/protein.faa > OUTPUT_FILENAME.faa

# Cluster protein seqeunces using MMseqs2

/home/jsrim/mmseqs/bin/mmseqs easy-cluster OUTPUT_FILENAME.faa result tmp --min-seq-id 0.5 --threads 55
```

3. Check _dcw_ gene presence in LAB genomes with custom script
4. Visulaization of _dcw_ genes present in the LAB genome with tools such as R [gggenomes](https://github.com/thackl/gggenomes), Python [pyGenomeViz](https://github.com/moshi4/pyGenomeViz)
