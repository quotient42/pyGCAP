#!/bin/bash

# 사용자로부터 TAXON 값을 입력 받음
read -p ">> Enter TAXON value: " TAXON
ROOT_DIR="./ncbi_dataset"
DATA_DIR="${ROOT_DIR}/data"
INPUT_DIR="${ROOT_DIR}/input"

# echo "<< Downloading genome data..."
# ./datasets download genome taxon "$TAXON" --reference --include protein,gff3,gbff --filename "${TAXON}.zip"

# echo "<< Unzipping files..."
# unzip "${TAXON}.zip" -d .

# mv $DATA_DIR $INPUT_DIR
# mkdir $DATA_DIR

# echo "<< Formatting data..."
# ./dataformat tsv genome \
#   --inputfile "${INPUT_DIR}/assembly_data_report.jsonl" \
#   --fields "accession,organism-name,assmstats-total-sequence-len,checkm-completeness,checkm-contamination,annotinfo-featcount-gene-total" \
#   > "${DATA_DIR}/assembly_report.tsv"

# echo "<< Concatenating protein files..."
# cat ${INPUT_DIR}/GC*/protein.faa > "${DATA_DIR}/${TAXON}.faa"

# echo "<< Clustering protein sequences..."
# mmseqs easy-cluster "${DATA_DIR}/${TAXON}.faa" "${DATA_DIR}/result" "${DATA_DIR}/tmp" --min-seq-id 0.5 --threads 55

echo "<< Organizing directories..."
mkdir ${ROOT_DIR}/output
mkdir ${ROOT_DIR}/output/genus
mkdir ${ROOT_DIR}/output/img
mkdir ${ROOT_DIR}/output/tsv
mkdir ${ROOT_DIR}/seqlib

# echo "<< Data organization complete."
