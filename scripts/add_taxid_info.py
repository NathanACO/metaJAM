import re
import pandas as pd
import numpy as np
import sys
import argparse
import time
#You need to install the taxadb package with "pip install taxadb" in order to import the 2 following modules:
from taxadb.taxid import TaxID
from taxadb.accessionid import AccessionID


# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Add taxonomy informations into blast tab.\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allows you to add taxonomic information to a blast, diamond or mmseqs2 output file. ')
parser.add_argument("-b", "--blast", help="The blast file in tabular format \t")
parser.add_argument("-o", "--out",help="The ouptut path (with the new filename) where to create the oufile")
parser.add_argument("-d", "--taxadb_sqlite_file", help="The directory where is located the sqlite database")
parser.add_argument("-blst", "--blast_type", help="The type of blast : support blast (blast), diamond (diamond) and mmseqs2 (mmseqs2) analysis")
parser.add_argument("-del", "--Output_delimitator",help="The desired output delimitator, by default : ';'  but you cand add tabular format with the option : 'tab'")



args = parser.parse_args()

"""
Ex usage python3 add_taxid_info.py -blst Mmseqs2  -b /beegfs/data/bguinet/these/Horizon_project_part/tBLASTn_VLPs/Matches_VLPs_prot_vs_NR_mmseqs2_all.m8 -d /beegfs/data/bguinet/taxadb/taxadb_new.sqlite -o /beegfs/data/bguinet/these/Horizon_project_part/tBLASTn_VLPs/Matches_VLPs_prot_vs_NR_mmseqs2_all_taxid.m8 
python3 python_test.py -blst Mmseqs2  -b tab_test -d /beegfs/data/bguinet/taxadb/taxadb_new.sqlite -o tab_test_taxid.m8 
"""

"""
Here are the steps to install the taxadb sqlite (from https://github.com/HadrienG/taxadb)

1) Install the package taxadb :
>>> pip3 install taxadb

2) Create the SQlite database:

(In order to create your own database, you first need to download the required data from the ncbi ftp)
>>> /beegfs/data/id/myconda/bin/taxadb download -o taxadb 

(If you work only on protein sequences, no need to download the all db (-d prot option only get protein informations)
>>> nohup /beegfs/data/id/myconda/bin/taxadb create -i /beegfs/data/id/taxadb --dbname taxadb.sqlite -d prot --fast &> nohup.out&

/!\ - Careful
There is a typo in the taxadb code in the taxid.py file, so you have to do a little nano on the folder and replace on the lineage_name function (line 107 lineages.append((rank, current_lineage_id)) by lineages.append((rank, current_lineage))
The taxid.py file shoud be located into your python lib environment  ex: /lib/python3.7/site-packages/taxadb/taxid.py

#The taxadb.sqlite file will by created in the directory where the command is ran.
"""

"""
Input exemple file:
query	target	pident	alnlen	qstart	qend	tstart	tend	evalue	bits
scaffold_1	ACT53096.1	33.3	102	1	295	158	259	1.19e-12	64
scaffold_2	ADR03158.1	28.9	183	3	531	98	271	2.516e-16	78
scaffold_3	BAV91453.1	26.1	168	9	486	63	227	2.294e-07	51

Output exemple file with tab delimitator:
	query	target	pident	alnlen	qstart	qend	tstart	tend	evalue	bits	Taxid	species	genus	subfamily	family	superfamily	parvorder	infraorder	suborder	order	superorder	no_rank	cla
ss	superclass	subphylum	phylum	kingdom	superkingdom	subgenus	subcohort	cohort	infraclass	subclass
0	scaffold_1	ACT53096.1	33.3	102	1	295	158	259	1.19e-12	64	9606	Homo_sapiens	Homo	Homininae	Hominidae	Hominoidea	Catarrhini	Simiiformes	Haplorrhini
	Primates	Euarchontoglires		cellular_organisms	Mammalia	Sarcopterygii	Craniata	Chordata	Metazoa	Eukaryota	NA	NA	NA	NA	NA
1	scaffold_2	ADR03158.1	28.9	183	3	531	98	271	2.516e-16	78	10090	Mus_musculus	Mus	Murinae	Muridae	NA	NA	NA	Myomorpha	Rodentia	Euarchontoglires	cellular_organisms	Mammalia	Sarcopterygii	Craniata	Chordata	Metazoa	Eukaryota	Mus	NA	NA	NA	NA
2	scaffold_3	BAV91453.1	26.1	168	9	486	63	227	2.294e-07	51	487618	Danio_margaritatus	Danio	NA	Cyprinidae	NA	NA	NA	Cyprinoidei	Cypriniformes	Cyp
riniphysae	cellular_organisms	Actinopteri	Actinopterygii	Craniata	Chordata	Metazoa	Eukaryota	NA	Ostariophysi	Otomorpha	Teleostei	Neopterygii
"""


#In order to get the processing time
start_time = time.time()

# Variable that stores the options given by the user
blast_file=args.blast
out_file=args.out
taxadb_sqlite_file=args.taxadb_sqlite_file
blast_type=args.blast_type

#blast_file="/beegfs/data/bguinet/these/Clustering/mmseqs2_search_analysis/Viralprot_vs_Viral_loci_result_all_match_and_cluster_without_duplicated.m8"
#taxadb_sqlite_file="/beegfs/data/bguinet/taxadb.sqlite"
#out_file="/beegfs/data/bguinet/these/Clustering/mmseqs2_search_analysis/Viralprot_vs_Viral_loci_result_all_match_and_cluster_without_duplicated_taxid.m8"

#Open the blast dataframe by keeping the column names 
blast=pd.read_csv(blast_file,sep="\t")

#remove all hits abox 0.01 evalue 

blast = blast.loc[blast['evalue'].lt(0.01)]

#print(blast)
#blast=blast.groupby('query').head(3) #in order to only have 3 estimations for taxid/queries to save time
print("\n")
print("#################################################")
print("Recovery process of taxonomic IDs in progress" )
print("#################################################")
print("\n")
# Create accession object
accession = AccessionID(dbtype='sqlite', dbname=taxadb_sqlite_file)

#Create a dataframe to store the Acc_number and the Taxid informations
blast_taxid=pd.DataFrame(columns=['Acc_number','Taxid'])


Acc_number_list=[]
if blast_type == 'blast' or blast_type == 'diamond' or blast_type == 'Blast' or blast_type == 'Diamond' :
	for i in blast.ssid.str.split(".").str[0]:    
	        Acc_number_list.append(i)
if blast_type == 'Mmseqs2' or blast_type == 'mmseqs2' or blast_type == 'mmseqs' or blast_type == 'Mmseqs':
	for i in blast.target.str.split(".").str[0]: 
	        Acc_number_list.append(i)
#Remove duplicate acc_number
Acc_number_list = list(dict.fromkeys(Acc_number_list))

#Taxadb has some problems processing refseq identifiers with the addition of a .number, so we will remove these numbers. 
Acc_number_list=[re.sub(r'\.\d+', '', i) for i in Acc_number_list]

#Function to return smaller lists of a larger list in order to question the database with several accession number instead of just one.
def chunk(items,n):
    for i in range(0,len(items),n):
        yield items[i:i+n]

number_of_sub_acc_number_lists=round(len(Acc_number_list)/990)+1
print(blast_taxid)

# Break the larger list into 10 length lists.
count_row = len(Acc_number_list)
filecount = 1

#For each sub-list of access numbers of size 10, retrieve the equivalent taxonomic IDs and stor it into taxids variable.
for sublist in chunk(Acc_number_list,number_of_sub_acc_number_lists):
	taxids = accession.taxid(sublist)
	#For each taxonomic ID in the variable, add in the blast_taxid dataframe the accession number in col1 and the corresponding taxonomic number in col2. 
	for tax in taxids:
		blast_taxid.loc[len(blast_taxid)]=[tax[0],tax[1]]
		#A script to add process progression bar
		filled_len = int(round(50 * filecount / float(count_row -1)))
		percents = round(100.0 * filecount / float(count_row ), 1)
		bar = '=' * filled_len + '-' * ((50) - filled_len)
		sys.stdout.write('[%s] %s%s Get Tax Id informations...%s\r' % (bar, percents, '%', ''))
		sys.stdout.flush()  # As suggested by Rom Ruben
		filecount += 1


# Add to the blast table the 'Taxid' column.
if blast_type == 'blast' or blast_type == 'diamond' or blast_type == 'Blast' or blast_type == 'Diamond' :
	blast['Taxid']=blast.ssid.str.split(".").str[0].map(dict(zip(blast_taxid.Acc_number,blast_taxid.Taxid)))
	blast['New_Acc_number'] =  [re.sub(r'[\n\r]*','', str(x)) for x in blast['ssid']]
if blast_type == 'Mmseqs2' or blast_type == 'mmseqs2' or blast_type == 'mmseqs' or blast_type == 'Mmseqs':
	blast['Taxid']=blast.target.str.split(".").str[0].map(dict(zip(blast_taxid.Acc_number,blast_taxid.Taxid)))
	blast['New_Acc_number'] =  [re.sub(r'[\n\r]*','', str(x)) for x in blast['target']]

print("Recovery process of taxonomic Ids done" )

print("\n")
print("########################################################")
print("Recovery process of taxonomic informations in progress" )
print("########################################################")
print("\n")

# Create taxid object
taxid = TaxID(dbtype='sqlite', dbname=taxadb_sqlite_file)

# Create a dictionary that will store the taxonomic informations corresponding to each taxonomic identifier 
Taxid_dictionnary = {}

count_row = blast_taxid.shape[0]
filecount = 1

# Pour chaque taxid dans le tableau, ajouter au dictionnaire d 
for tax_ID in blast_taxid['Taxid']:
	# Query the SQlite database and store the result into the dictionnary
	Taxid_dictionnary[tax_ID] = taxid.lineage_name(tax_ID, ranks=True,reverse=True)
	#A script to add process progression bar
	filled_len = int(round(50 * filecount / float(count_row -1)))
	percents = round(100.0 * filecount / float(count_row ), 1)
	bar = '=' * filled_len + '-' * ((50) - filled_len)
	sys.stdout.write('[%s] %s%s Get taxonomic informations...%s\r' % (bar, percents, '%', ''))
	sys.stdout.flush()  # As suggested by Rom Ruben
	filecount += 1

print("Recovery process of taxonomic informations done" )


#Add taxonomic informations into the blast dataframe
blast=blast.merge(pd.DataFrame({k: dict(v) for k, v in Taxid_dictionnary.items()}).T, left_on='Taxid', right_index=True,how='outer')
#Remove NA values
blast.fillna('NA', inplace=True)
#Remove all spaces en replace them by _ - double check
blast.columns = blast.columns.str.replace(' ', '_')
blast=blast.replace(' ', '_', regex=True)
#Remove the column 'New_Acc_number'
blast=blast.drop(columns=['New_Acc_number'])

"""
#################
#ADDD genomic structure informations 

Virus_genomic_structure_tab=pd.read_csv("/beegfs/data/bguinet/these/NCBI_nucleotides_viruses/All_virus_families_informations.csv",sep=";",header=0,encoding='latin-1')
###This table comes from the  Virus Metadata Repository: version November 27, 2019; MSL34
Virus_genomic_structure_tab=Virus_genomic_structure_tab[['Family','Genome composition']]
Virus_genomic_structure_tab = Virus_genomic_structure_tab.drop_duplicates()

Virus_genomic_structure_tab.columns = ['family', 'genomic_structure']
Virus_genomic_structure_tab=Virus_genomic_structure_tab.dropna()


blast =pd.merge(blast, Virus_genomic_structure_tab, left_on='family', right_on='family',how='outer')
#blast = blast[blast['Clustername'].notna()]
#Carrefull here in the customised code only representative of the result of this paper
blast.loc[blast['species'] == "Leptopilina_boulardi_filamentous_virus", ['genomic_structure']] = 'dsDNA'
blast.loc[blast['species'] == "Leptopilina_boulardi_filamentous_virus", ['family']] = 'LbFV-like family'
blast.loc[blast['target'].str.contains("QKN22",na=False),['genomic_structure']] = 'dsDNA' 
blast.loc[blast['target'].str.contains("QKN22",na=False),['species']] = "Drosophila-associated filamentous virus"
blast.loc[blast['target'].str.contains("QKN22",na=False),['family']] = 'LbFV-like family' 


#Remove NA 
blast = blast[blast['query'].notna()]
"""

blast.to_csv(out_file,sep=";",index=False)

#Merge all information with the revious blast dataframe 
"""
blast2=pd.read_csv(blast_file,header=0,sep=";")

#blast.drop(['Unnamed:_0', 'Clustername', 'query', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tlen'], axis = 1, inplace = True) 

blast.drop(['query', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits'], axis = 1, inplace = True) 
blast3=pd.merge(blast2,blast,on='target', how='outer')
blast3=blast3.drop_duplicates(subset=['query', 'target'], keep='first')



print("Number of queries with a family assigned : ", len(blast3[blast3['family'].notna()]['query'].unique()))
print("Number of queries with a species assigned : ", len(blast3[blast3['species'].notna()]['query'].unique()))

if args.Output_delimitator=="tab":
	blast3.to_csv(out_file, sep="\t",index=False)
else:
	blast3.to_csv(out_file, sep=";",index=False)
"""


print("\n")
print("Process done")
print("Output written at: ",out_file)
print("--- %s seconds ---" % (time.time() - start_time))

