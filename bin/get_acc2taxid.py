#!/usr/bin/env python3
from Bio import Entrez
import sys
import time
arg1 = sys.argv[1]

Entrez.email = arg1+"@gmail.com"

MAX_RETRY = 6
SLEEP = 0.3

def get_taxid(species):
    for _ in range(MAX_RETRY):
        try:
            handle = Entrez.esearch(db="taxonomy", term=species)
            record = Entrez.read(handle)
            if record["IdList"]:
                return record["IdList"][0]
        except Exception:
            pass
        time.sleep(SLEEP)
    return "NA"

#initiate a dictionary with speices and taxid
dict_species_taxid = {}

with open(arg1) as f:
    for line in f:
        species = " ".join(line.strip().split())

        #if the species is not recorded
        if species not in dict_species_taxid:
            taxid = get_taxid(species)
            dict_species_taxid[species] = taxid
        else:
            taxid = dict_species_taxid[species]

        print(f"{species}\t{taxid}")
        sys.stdout.flush()
