import re
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from time import sleep
from pathlib import Path
import argparse
import openpyxl as op
import os

# Configuration
Entrez.email = 'jp.uesc17@gmail.com'
Entrez.api_key = '7533e87f7eb6af740cc79274b033e9427908'

# Load reference files
sequences_df = pd.read_csv('sequences.tsv', sep='\t')
sequences_df.drop_duplicates(subset="Accession", inplace=True)

local_virus_file = 'ICTV_Master_Species_List_2022_MSL38.v2.xlsx'
ictv_df = pd.read_excel(local_virus_file, sheet_name='MSL')
family_col = ictv_df['Family']
genome_col = ictv_df['Genome Composition']

# Global dictionaries
data_map = {}
genome_map = {}

def get_ncbi_taxonomy(taxon):
    try:
        # If taxon is not a numeric ID, search for it first
        if not re.match(r'^\d+$', taxon):
            handle = Entrez.esearch(db='taxonomy', term=f'"{taxon}"', rettype='gb', retmode='text')
            record = Entrez.read(handle, validate=False)
            handle.close()
            sleep(0.5)
            
            if record['IdList']:
                tax_id = record['IdList'][0]
            else:
                tax_id = input(f"What is the ID code for {taxon}?\n")
        else:
            tax_id = taxon

        # Fetch taxonomy details using ID
        handle2 = Entrez.efetch(db='taxonomy', id=tax_id, retmode='xml')
        record2 = Entrez.read(handle2, validate=False)
        handle2.close()
        sleep(0.5)

        # Extract Family from lineage
        for tax_element in reversed(record2[0]['LineageEx']):
            if tax_element['Rank'] == 'family':
                return tax_element['ScientificName']
        return 'Family not found'
    
    except Exception as e:
        print(f"Error fetching taxonomy for {taxon}: {e}")
        sleep(60)
        return get_ncbi_taxonomy(taxon)

def extract_protein_code(protein):
    matches = re.findall(r'(\w+)', str(protein))
    return matches[0] if matches else 'Unknown'

def update_local_database():
    global data_map, genome_map
    local_db_file = 'local_database.csv' # Renamed from banco_dados_local.csv
    if Path(local_db_file).exists():
        db_df = pd.read_csv(local_db_file)
        for _, row in db_df.iterrows():
            data_map[row['Organism_Name']] = row['Family']
            genome_map[row['Family']] = row['Genome']

def save_local_database():
    db_df = pd.DataFrame([
        {'Organism_Name': esp, 'Family': fam, 'Genome': genome_map.get(fam, 'Unknown')} 
        for esp, fam in data_map.items()
    ])
    db_df.to_csv('local_database.csv', index=False)
    print("Local database saved as local_database.csv")

def process_family(file_path, accession):
    global data_map, genome_map
    table = pd.read_table(file_path, header=None, sep="\t")
    protein_codes = [extract_protein_code(p) for p in table[2]] 

    family_list, genome_list = [], []
    
    for code in protein_codes:
        if code in data_map:
            family_list.append(data_map[code])
        elif code in sequences_df["Accession"].values:
            molecule_type = sequences_df.loc[sequences_df["Accession"] == code, "Molecule_type"].values[0]
            data_map[code] = molecule_type
            family_list.append(molecule_type)
        else:
            print(f"Fetching: {code}")
            sleep(4)
            data_map[code] = get_ncbi_taxonomy(code)
            family_list.append(data_map[code])

    for family in family_list:
        if family in genome_map:
            genome_list.append(genome_map[family])
        else:
            if family in family_col.values:
                index = family_col[family_col == family].index[0]
                genome_map[family] = genome_col[index]
            else:
                genome_map[family] = family
            genome_list.append(genome_map[family])

    records = [
        SeqRecord(
            Seq(seq), 
            id=protein_codes[i], 
            description=f"Family: {family_list[i]}, Genome: {genome_list[i]}"
        ) 
        for i, seq in enumerate(table[2])
    ]
    
    # 📁 Output path
    output_folder = f"EVEs_of_{accession}" # Renamed folder structure slightly
    os.makedirs(output_folder, exist_ok=True)

    annotated_fasta = os.path.join(output_folder, f"{accession}_blastx_annotated.fasta")
    fasta_no_rt = os.path.join(output_folder, f"{accession}_no_RT.fasta")

    SeqIO.write(records, annotated_fasta, "fasta")
    print(f"FASTA file saved as {annotated_fasta}")

    records_no_rt = [r for r in records if 'RT' not in r.description.upper()]
    SeqIO.write(records_no_rt, fasta_no_rt, "fasta")
    print(f"File without RT saved as {fasta_no_rt}")
    
def filter_blast_results(lines):
    VIRAL_KEYWORDS = [
        "virus", "viral", "virion", "capsid", "coat", "envelope", "spike",
        "helicase", "polymerase", "transcriptase", "replicase", "viroporin",
        "attachment", "fusion", "ns1", "ns2", "ns3", "nsp", "vpg", "reverse transcriptase",
        "polyprotein", "hypothetical protein", "rdrp"
    ]
    # Note: These names are scientific, usually kept in Latin/English, so no change needed
    EXCLUDE_ORGANISMS = [
        "acanthamoeba castellanii medusavirus",
        "mollivirus sibericum",
        "pacmanvirus a23",
        "kaumoebavirus"
    ]

    input_db = [l.strip().split('\t') for l in lines]
    print(f"Total entries: {len(input_db)}")
    
    # Logic to keep the best hit (lowest E-value/Higher Score) per query
    output_db = [input_db[0]]

    for i in input_db[1:]:
        if i[0] == output_db[-1][0]:
            # Comparing E-value (index 5) and Bitscore/Identity (index 6)
            if float(i[5]) < float(output_db[-1][5]) or (float(i[5]) == float(output_db[-1][5]) and float(i[6]) > float(output_db[-1][6])):
                output_db[-1] = i
        else:
            output_db.append(i)

    filtered_db = [output_db[0]]

    for row in output_db[1:]:
        col_split = [field.strip() for field in row[4].split("|")]
        
        # Ensure the column exists before accessing
        if len(col_split) > 5:
            protein_name = col_split[1].lower()
            organism = col_split[3].lower()
            genome_type = col_split[5].lower()

            if any(g in genome_type for g in ["dsdna", "dsdna-rt", "ssrna-rt"]):
                continue
            else:
                filtered_db.append(row)
        else:
            # Fallback if split format is unexpected
            filtered_db.append(row)

    print(f"Total after filtering: {len(filtered_db)}")
    return filtered_db

def save_tsv(data, name):
    tsv_name = name.replace('.xlsx', '.tsv') if name.endswith('.xlsx') else name
    df = pd.DataFrame(data[1:], columns=data[0])
    df.to_csv(tsv_name, sep='\t', index=False)
    print(f"TSV file saved as {tsv_name}")

def main():
    parser = argparse.ArgumentParser(description="Processes input files and updates the database.")
    parser.add_argument("diamond_out", help="Input file (DIAMOND output) to process families")
    parser.add_argument("input_file", help="Input file for extraction (e.g., full blastx)")
    parser.add_argument("output_tsv", help="Name of the output TSV file")
    args = parser.parse_args()

    # 🔍 Extract accession from DIAMOND filename
    accession = os.path.basename(args.diamond_out).split("_")[0]
    output_folder = f"EVEs_of_{accession}"
    os.makedirs(output_folder, exist_ok=True)

    update_local_database()

    # 📁 Automatically adjusts annotated FASTA output path
    process_family(args.diamond_out, accession)

    save_local_database()

    with open(args.input_file) as f:
        lines = f.readlines()

    filtered_db = filter_blast_results(lines)

    # Saves the TSV inside the correct folder
    tsv_name = os.path.join(output_folder, args.output_tsv)
    save_tsv(filtered_db, tsv_name)

if __name__ == "__main__":
    main()
