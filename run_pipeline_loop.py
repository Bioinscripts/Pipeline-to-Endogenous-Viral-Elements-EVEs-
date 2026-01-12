import subprocess
from pathlib import Path
import sys

# File containing genome accessions
genome_list_file = "genomes.txt"

# Fixed database paths
# IMPORTANT: Update these paths to your actual English paths if they changed
protein_db = "/home/joaopedro/Downloads/sequences_1.fasta"
nucleotide_db = "/home/joaopedro/work_local/Laiana/seq/data_cluster_laiana.fasta"
threads = "8"

# Create log directory
Path("logs").mkdir(exist_ok=True)

# Read accessions
if not Path(genome_list_file).exists():
    sys.exit(f"Error: {genome_list_file} not found.")

with open(genome_list_file) as f:
    accessions = [line.strip() for line in f if line.strip()]

for acc in accessions:
    print(f"🔁 Processing {acc}...")

    log_file_path = f"logs/{acc}.log"
    
    # We use 'eve_pipeline.py' assuming you renamed the main script
    with open(log_file_path, "w") as log_file:
        process = subprocess.run([
            "python", "eve_pipeline.py",
            "--c", acc,
            "--p", protein_db,
            "--n", nucleotide_db,
            "--t", threads,
            "--download" # Assuming you want to download based on the loop logic
        ], stdout=log_file, stderr=subprocess.STDOUT)

    if process.returncode == 0:
        print(f"[✓] {acc} finished successfully.")
    else:
        print(f"[X] Error processing {acc}. Check log: {log_file_path}")
