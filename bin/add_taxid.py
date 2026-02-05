#!/usr/bin/env python3
import sys

# input_contig = sys.argv[1]
# species_taxid_file = sys.argv[2]

input_contig = "mammals_birds_fish.contig"
species_taxid_file = "species_taxid"

# read species_taxid into dict: species -> taxid
species_taxid = {}
with open(species_taxid_file) as f:
    for line in f:
        if not line.strip():
            continue
        species, taxid = line.rstrip().split('\t', 1)
        species_taxid[species] = taxid

# read input_contig, map second column to taxid
with open(input_contig) as f:
    for line in f:
        contig, *species1 = line.strip().split()
        species = " ".join(species1[0:2])
        taxid = species_taxid[species]
        if taxid:
            print(contig, taxid)
