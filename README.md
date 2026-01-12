# EVE Discovery Pipeline

## Overview
This repository contains a bioinformatics pipeline designed to identify, annotate, and analyze **Endogenous Viral Elements (EVEs)** derived from non-retroviral RNA viruses within eukaryotic genomes (with a focus on fungi, though fully adaptable).

The pipeline automates genome acquisition, ORF prediction, homology-based filtering, viral annotation, redundancy removal, and taxonomic classification using NCBI and ICTV resources.

## 🚀 Features
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
3. **DIAMOND**
4. **CD-HIT**
5. **ORFipy**
