import time
import csv
import io
import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML


def run_blast(batch_seqs, program="blastn", db="nt", num_alignments=20, max_retries=5):
    """
    Run an online BLAST search at NCBI for a batch of sequences,
    using progressive backoff between retries.

    Retry schedule:
        1st attempt: immediate
        2nd attempt: 5 minutes
        3rd attempt: 15 minutes
        4th attempt: 30 minutes
        5th attempt: 1 hour

    Parameters:
        batch_seqs (list): List of tuples (ID, sequence)
        program (str): BLAST program (blastn, blastx, etc.)
        db (str): NCBI database (nt, nr, etc.)
        num_alignments (int): Maximum number of hits per sequence
        max_retries (int): Maximum number of retries
    """
    backoff_times = [0, 300, 900, 1800, 3600]  # seconds

    fasta_str = "\n".join([f">{header}\n{seq}" for header, seq in batch_seqs])

    for attempt in range(max_retries + 1):
        try:
            if attempt > 0:
                wait_time = backoff_times[min(attempt, len(backoff_times) - 1)]
                print(f"[!] Waiting {wait_time // 60} minutes before retrying...")
                time.sleep(wait_time)

            print(f"[{program.upper()}] Attempt {attempt + 1} of {max_retries + 1}")
            result_handle = NCBIWWW.qblast(
                program,
                db,
                fasta_str,
                hitlist_size=num_alignments,
                format_type="XML"
            )
            return result_handle.read()

        except Exception as e:
            print(f"[X] Error on attempt {attempt + 1}: {e}")
            if attempt == max_retries:
                print("[X] Maximum number of retries exceeded. Aborting BLAST.")
                return None


def highlight_mismatches(query_seq, match_seq, subject_seq):
    """Return query sequence highlighting mismatches and gaps."""
    highlighted = ""
    for q, m, s in zip(query_seq, match_seq, subject_seq):
        if q == "-":
            highlighted += "-"
        elif q != s:
            highlighted += q.lower()
        else:
            highlighted += q
    return highlighted


def parse_blast_results(xml_result):
    """Extract relevant information from BLAST XML results."""
    results = []
    alignments_data = []

    xml_handle = io.StringIO(xml_result)
    blast_records = NCBIXML.parse(xml_handle)

    for record in blast_records:
        query_id = record.query
        for alignment in record.alignments:
            subject_id = alignment.hit_def
            for hsp in alignment.hsps:
                query_coverage = (hsp.align_length / record.query_length) * 100
                highlighted_query = highlight_mismatches(
                    hsp.query, hsp.match, hsp.sbjct
                )

                results.append([
                    query_id,
                    subject_id,
                    alignment.length,
                    hsp.identities,
                    hsp.align_length,
                    record.query_length,
                    query_coverage,
                    hsp.query_start,
                    hsp.query_end,
                    hsp.sbjct_start,
                    hsp.sbjct_end,
                    hsp.expect
                ])

                alignments_data.append([
                    query_id,
                    subject_id,
                    hsp.query_start,
                    hsp.query_end,
                    hsp.sbjct_start,
                    hsp.sbjct_end,
                    hsp.expect,
                    highlighted_query,
                    hsp.match,
                    hsp.sbjct
                ])

    return results, alignments_data


def save_results_to_tsv(results, filename):
    """Save BLAST summary results to a TSV file."""
    headers = [
        "Query ID", "Subject ID", "Subject Length", "Identities",
        "Alignment Length", "Query Length", "Query Coverage",
        "Query Start", "Query End", "Subject Start", "Subject End", "E-value"
    ]
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(results)


def save_alignments_to_tsv(alignments, filename):
    """Save detailed alignments (with graphical representation) to TSV."""
    headers = [
        "Query ID", "Subject ID",
        "Query Start", "Query End",
        "Subject Start", "Subject End",
        "E-value",
        "Query Sequence", "Alignment", "Subject Sequence"
    ]
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(alignments)


def batch_blast(fasta_file, output_dir=None, batch_size=5):
    accession = os.path.basename(fasta_file).split("_clusters")[0]

    if output_dir is None:
        output_dir = f"EVEs_from_{accession}"

    os.makedirs(output_dir, exist_ok=True)

    sequences = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(fasta_file, "fasta")]

    all_blastn_results = []
    all_blastx_results = []
    all_blastn_alignments = []
    all_blastx_alignments = []
    analyzed_sequences = set()

    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]

        print(f"Processing batch {i // batch_size + 1} - BLASTn")
        xml_result = run_blast(batch, "blastn")
        if xml_result:
            results, alignments = parse_blast_results(xml_result)
            all_blastn_results.extend(results)
            all_blastn_alignments.extend(alignments)
            analyzed_sequences.update(r[0] for r in results)

        time.sleep(5)

        print(f"Processing batch {i // batch_size + 1} - BLASTx")
        xml_result = run_blast(batch, "blastx", db="nr")
        if xml_result:
            results, alignments = parse_blast_results(xml_result)
            all_blastx_results.extend(results)
            all_blastx_alignments.extend(alignments)
            analyzed_sequences.update(r[0] for r in results)

        time.sleep(5)

    save_results_to_tsv(
        all_blastn_results, os.path.join(output_dir, "all_blastn_results.tsv")
    )
    save_alignments_to_tsv(
        all_blastn_alignments, os.path.join(output_dir, "all_blastn_alignments.tsv")
    )
    save_results_to_tsv(
        all_blastx_results, os.path.join(output_dir, "all_blastx_results.tsv")
    )
    save_alignments_to_tsv(
        all_blastx_alignments, os.path.join(output_dir, "all_blastx_alignments.tsv")
    )

    all_ids = {rec.id for rec in SeqIO.parse(fasta_file, "fasta")}
    unanalyzed_ids = all_ids - analyzed_sequences

    with open(os.path.join(output_dir, "unanalyzed_sequences.fasta"), "w") as out_fa:
        for rec in SeqIO.parse(fasta_file, "fasta"):
            if rec.id in unanalyzed_ids:
                SeqIO.write(rec, out_fa, "fasta")

    print("Sequence alignment finished successfully.")

