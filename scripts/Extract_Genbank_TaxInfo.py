import re
import pandas as pd
import numpy as np
import sys
import argparse
import time 
from Bio import Entrez
from tqdm import tqdm

api_key="ce154668962f6a9d4763cab6df03b7add608"
Entrez.email = "benjamin.guinet95@gmail.com"  # Replace with your email
Entrez.api_key = api_key

print('-----------------------------------------------------------------------------------\n')
print('                        Add taxonomy information into blast tab.\n')
print('-----------------------------------------------------------------------------------\n')

parser = argparse.ArgumentParser(description='Add taxonomic information to a BLAST output file.')
parser.add_argument("-b", "--blast", help="The BLAST file in tabular format")
parser.add_argument("-o", "--out", help="The output file path")
parser.add_argument("-blst", "--blast_type", help="The type of BLAST (blast, diamond, mmseqs2)")
args = parser.parse_args()

start_time = time.time()

blast_file = args.blast
out_file = args.out
blast_type = args.blast_type


""" 
eg usage python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Extract_Genbank_TaxInfo.py -b /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Herpes_capture/New_capture/ALL_P34656_548689_vs_NT.m8 -o /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/Herpes_capture/New_capture/ALL_P34656_548689_vs_NT_taxo.m8 -blst mmseqs2
""" 
df_blast = pd.read_csv(blast_file, sep="\t")
df_blast = df_blast[df_blast['evalue'] < 0.01]  # Filter hits based on e-value

print("\n#################################################")
print("Recovery process of taxonomic IDs and taxo info in progress")
print("#################################################\n")

accession_list = df_blast['target'].unique()

def get_taxid_from_accessions(accessions):
    try:
        handle = Entrez.esearch(db="nucleotide", term=" OR ".join(accessions), retmax=len(accessions))
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            return {acc: None for acc in accessions}
        genbank_ids = record["IdList"]
        handle = Entrez.efetch(db="nucleotide", id=",".join(genbank_ids), rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        results = {rec.get("GBSeq_primary-accession", None): rec.get("GBSeq_taxonomy", None) for rec in records}
        return results
    except Exception as e:
        print("Error:", e)
        return {acc: None for acc in accessions}

batch_size = 500
batches = [accession_list[i:i + batch_size] for i in range(0, len(accession_list), batch_size)]
results_dict = {}

for batch in tqdm(batches, desc="Processing Batches (100)"):
    batch_results = get_taxid_from_accessions(batch)
    results_dict.update(batch_results)

df_blast.to_csv(out_file, sep="\t", index=False)

failed_accessions = [acc for acc, tax in results_dict.items() if tax is None]

if failed_accessions:
    batch_size = 10
    failed_batches = [failed_accessions[i:i + batch_size] for i in range(0, len(failed_accessions), batch_size)]
    for batch in tqdm(failed_batches, desc="Retrying Failed (10)"):
        results_dict.update(get_taxid_from_accessions(batch))

results_df = pd.DataFrame.from_dict(results_dict, orient="index", columns=["Taxonomy"]).reset_index()
results_df.rename(columns={"index": "Accession"}, inplace=True)

df_blast = df_blast.merge(results_df, left_on="target", right_on="Accession", how="left")

df_blast.to_csv(out_file, sep="\t", index=False)
