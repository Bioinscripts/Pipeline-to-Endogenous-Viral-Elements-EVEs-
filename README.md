# EVE Discovery Pipeline

## Introduction
This repository contains a bioinformatics pipeline designed to identify, annotate, and analyze **Endogenous Viral Elements (EVEs)** derived from non-retroviral RNA viruses within eukaryotic genomes (with a focus on fungi, though fully adaptable).

The pipeline automates genome acquisition, ORF prediction, homology-based filtering, viral annotation, redundancy removal, and taxonomic classification using NCBI and ICTV resources.

## Features
- **Automated Genome Downloading**: Retrieves genomes directly from NCBI using accession numbers.
- **ORF Prediction**: Detects Open Reading Frames using `orfipy`.
- **Host Filtering**: Removes host-like sequences via local BLASTn.
- **Viral Homology Search**: Performs fast protein alignments using DIAMOND BLASTx.
- **Taxonomic Assignment**: Retrieves lineage data from NCBI Taxonomy and maps results to ICTV classifications.
- **Redundancy Removal**: Clusters sequences using CD-HIT-EST to generate a non-redundant dataset.
- **Online Validation**: Performs batch BLASTn/BLASTx searches against NCBI remote databases for final confirmation.

## Prerequisites

Ensure the following tools are installed and available in your system path:

1. **Python 3.8+**
2. **NCBI BLAST+** (`makeblastdb`, `blastn`)
3. **DIAMOND** (`blastx`)
4. **CD-HIT**
5. **ORFipy**
6. **ENTREZ**
7. **Biopython**


## Install Databases

### RefSeq Viral release
Podemos baixar as proteínas virais para o alinhamento a nível de blastX com as ORFs preditas em cada genoma através do site: (`https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/`).
Utiliza-se o banco de dados proteico (`https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz`).

<img width="502" height="320" alt="image" src="https://github.com/user-attachments/assets/19cd1269-fb19-46df-9354-171194cd77c8" />

\**It is important to note that additional databases may be utilized according to the specific objectives of each study.**\

### Taxonomy
  
The taxonomic data were obtained from the ICTV table released in 2023 (`ICTV_Master_Species_List_2022_MSL38.v2.xlsx`). For newly described species or those not yet cataloged, taxonomy was inferred in conjunction with the NCBI database, using protein identifiers and their corresponding genetic material. Owing to the recurrent presence of certain proteins, a reference database was constructed from previous results to facilitate the initial identification of genetic material. This database was subsequently employed to advance the pipeline for less common sequences that required further investigation via ENTREZ.

### 
