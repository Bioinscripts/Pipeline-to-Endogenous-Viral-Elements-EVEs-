import os
from Bio import SeqIO

def load_ids(*files):
    ids = set()
    for file_path in files:
        if not os.path.exists(file_path):
            continue
        with open(file_path) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split("\t")
                    ids.add(parts[0])
    return ids

def extract_sequences(input_fasta, target_ids, output_fasta):
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in target_ids:
                SeqIO.write(record, out, "fasta")

def remove_redundancy_cd_hit(input_fasta, output_fasta, identity=0.9, coverage=0.9):
    cmd = (
        f"cd-hit-est -i {input_fasta} -o {output_fasta} -c {identity} -aS {coverage} -M 40000 -T 4"
    )
    os.system(cmd)

def cd_hit_filter_pipeline(accession, input_folder=None):
    # 📁 Define folder based on accession if not provided
    if input_folder is None:
        input_folder = f"EVEs_of_{accession}" # Matched the folder name change in the previous script

    blastn_path = os.path.join(input_folder, f"{accession}_orfs_blastn.tsv")
    blastx_path = os.path.join(input_folder, f"{accession}_blastx_filtered.tsv") # Renamed from _filtrado
    viral_fasta = os.path.join(input_folder, f"{accession}_viral_hits.fasta")
    all_orfs_fasta = f"{accession}_orfs/{accession}_orfs.fasta"

    ids = load_ids(blastn_path, blastx_path)

    if os.path.exists(viral_fasta):
        for record in SeqIO.parse(viral_fasta, "fasta"):
            ids.add(record.id)

    combined_fasta = os.path.join(input_folder, f"{accession}_combined_filtered.fasta")
    extract_sequences(all_orfs_fasta, ids, combined_fasta)

    cluster_out = os.path.join(input_folder, f"{accession}_clusters_nr.fasta")
    remove_redundancy_cd_hit(combined_fasta, cluster_out)

    print(f"[✓] Combined and non-redundant sequences saved at: {cluster_out}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--accession", required=True, help="Input file prefix (e.g., GCF_xxxx)")
    parser.add_argument("--input_folder", default="results", help="Folder containing input files")
    args = parser.parse_args()
    
    cd_hit_filter_pipeline(args.accession, args.input_folder)
