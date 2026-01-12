import requests
import subprocess
import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

# Updated imports to match the new English filenames and function names
from blast_processing import (
    filter_blast_results,        # formerly extracao
    save_tsv,
    process_family,              # formerly processar_familia
    update_local_database,       # formerly atualizar_banco_dados_local
    save_local_database          # formerly salvar_banco_dados_local
)
from remove_duplicates import cd_hit_filter_pipeline # formerly pipeline_filtro_cd_hit
from online_blast import batch_blast # Ensure the file is named online_blast.py


# ============================
# PIPELINE STEPS
# ============================

def download_genome_ncbi(accession, file_type="fna"):
    """
    Download a genome file from the NCBI FTP based on accession
    (e.g. GCF_000240135.3).

    Parameters:
        accession (str): Genome accession
        file_type (str): File type (fna, gbff, gff, etc.)
    """
    try:
        prefix, rest = accession.split("_")
        number, version = rest.split(".")

        blocks = [number[i:i + 3] for i in range(0, len(number), 3)]
        path = "/".join(blocks)

        asm_number = int(number)
        base_name = f"{accession}_ASM{asm_number}v{version}"
        file_name = f"{base_name}_genomic.{file_type}.gz"

        url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{path}/{accession}/{file_name}"
        print(f"Downloading from: {url}")

        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(file_name, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"{file_name} downloaded successfully.")
        else:
            print(f"Error {response.status_code}: file not found or access denied.")

    except Exception as e:
        print(f"Error downloading {accession}: {e}")


def predict_orfs(genome_fasta, accession):
    print(f"Reading input genome: {genome_fasta}")

    combined_fasta = f"{accession}_with_complement.fasta"

    with open(combined_fasta, "w") as out:
        for record in SeqIO.parse(genome_fasta, "fasta"):
            SeqIO.write(record, out, "fasta")

            comp_record = record[:]
            comp_record.seq = record.seq.complement()
            comp_record.id = record.id + "_complement"
            comp_record.description = "complement strand (same 5' to 3' direction)"
            SeqIO.write(comp_record, out, "fasta")

    print("Original and complement strands saved.")

    # Note: orfipy command remains the same as it's a CLI tool
    cmd = (
        f"orfipy {combined_fasta} "
        f"--out {accession}_orfs "
        f"--min 100 --max 6000 "
        f"--dna {accession}_orfs.fasta "
        f"--strand b --ignore-case --include-stop"
    )
    os.system(cmd)

    print("ORF prediction completed.")


def sequence_alignment(
    protein_db,
    nucleotide_db,
    orf_fasta,
    accession,
    threads=4,
    output_dir="results"
):
    os.makedirs(output_dir, exist_ok=True)

    # ---------- BLASTn ----------
    blastn_out = os.path.join(output_dir, f"{accession}_orfs_blastn.tsv")

    # Check if blastdb files exist (checking common extensions)
    if not all(os.path.exists(f"{nucleotide_db}.{ext}") for ext in ("nhr", "nin", "nsq")):
        print(f"[!] Nucleotide DB not indexed. Running makeblastdb...")
        os.system(f"makeblastdb -in {nucleotide_db} -dbtype nucl")

    os.system(
        f"blastn -query {orf_fasta} -db {nucleotide_db} -out {blastn_out} "
        f"-num_threads {threads} "
        f"-outfmt '6 qseqid sseqid pident length mismatch gapopen "
        f"qstart qend sstart send evalue bitscore qcovs slen' "
        f"-max_target_seqs 1"
    )

    aligned_ids = set()
    # Safely handle if blastn_out wasn't created (e.g., no hits or error)
    if os.path.exists(blastn_out):
        with open(blastn_out) as f:
            for line in f:
                # Check E-value column (index 10)
                parts = line.split("\t")
                if len(parts) > 10 and float(parts[10]) < 1e-5:
                    aligned_ids.add(parts[0])

    filtered_fasta = os.path.join(output_dir, f"{accession}_filtered.fasta")
    with open(orf_fasta) as inp, open(filtered_fasta, "w") as out:
        for record in SeqIO.parse(inp, "fasta"):
            if record.id not in aligned_ids:
                SeqIO.write(record, out, "fasta")

    print(f"[✓] {len(aligned_ids)} ORFs removed after host-genome alignment.")

    # ---------- DIAMOND blastx ----------
    blastx_out = os.path.join(output_dir, f"{accession}_orfs_diamond.tsv")

    # Creating Diamond DB (if needed, though makedb overwrites usually)
    subprocess.run(
        ["diamond", "makedb", "--in", protein_db, "-d", protein_db],
        check=True
    )

    subprocess.run([
        "diamond", "blastx",
        "--query", filtered_fasta,
        "--db", protein_db,
        "--out", blastx_out,
        "--outfmt", "6",
        "qseqid", "qlen", "sseqid", "slen", "stitle",
        "evalue", "bitscore", "qcovhsp", "qtitle",
        "--max-target-seqs", "5",
        "--evalue", "1e-10",
        "--threads", str(threads)
    ], check=True)

    # Process results using functions from blast_processing.py
    update_local_database()
    process_family(blastx_out, accession)
    save_local_database()

    with open(blastx_out) as f:
        lines = f.readlines()

    filtered_db = filter_blast_results(lines)
    save_tsv(
        filtered_db,
        os.path.join(output_dir, f"{accession}_blastx_filtered.tsv")
    )


def main():
    parser = argparse.ArgumentParser(
        description="Pipeline for EVE analysis in fungal genomes"
    )

    parser.add_argument("--download", action="store_true", help="Download genome from NCBI")
    parser.add_argument("--c", dest="accession", help="Genome accession (e.g., GCF_...)")
    parser.add_argument("--g", dest="genome_file", help="Local genome FASTA file")
    parser.add_argument("--p", dest="protein_db", required=True, help="Protein database path")
    parser.add_argument("--n", dest="nucleotide_db", required=True, help="Nucleotide database path")
    parser.add_argument("--t", dest="threads", type=int, default=4, help="Number of threads")

    args = parser.parse_args()

    if args.download and not args.accession:
        sys.exit("Error: --c (accession) is required when using --download")

    if not args.download and not args.genome_file:
        sys.exit("Error: --g (genome file) is required when not using --download")

    accession = args.accession
    
    # If downloading, genome name is derived from accession
    if args.download:
        download_genome_ncbi(accession)
        # Note: download_genome_ncbi saves as .gz, you might need to unzip it here or update logic
        # For this script, assuming the user unzips or the download function saves as .fasta
        # Since the original code assumed "accession.fasta", we adapt:
        genome_path = f"{accession}.fasta" 
        # WARNING: The download function saves as .gz. You may need to add gunzip logic here.
    else:
        genome_path = args.genome_file
        if not accession:
            accession = os.path.splitext(os.path.basename(genome_path))[0]

    output_dir = f"EVEs_of_{accession}" # Renamed to match other scripts
    os.makedirs(output_dir, exist_ok=True)

    predict_orfs(genome_path, accession)

    orf_fasta_path = f"{accession}_orfs/{accession}_orfs.fasta"

    sequence_alignment(
        args.protein_db,
        args.nucleotide_db,
        orf_fasta_path,
        accession,
        args.threads,
        output_dir
    )

    cd_hit_filter_pipeline(accession, output_dir)

    clustered_fasta = os.path.join(
        output_dir, f"{accession}_clusters_nr.fasta"
    )

    batch_blast(clustered_fasta, output_dir=output_dir)


if __name__ == "__main__":
    main()
