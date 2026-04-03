#!/bin/bash -l

input_contig=$1

cut -f 2,3 -d' ' $input_contig | sort | uniq > species

./get_acc2taxid.py species > species_taxid

./add_taxid.py $input_contig species_taxid > $input_contig.acc2taxid

