import pandas as pd
import argparse
from subprocess import call
import os

parser = argparse.ArgumentParser()
parser.add_argument("fasta", nargs='?', help="Input dataset fasta file")
parser.add_argument('--num_threads', help='Number of threads for blastn and seqkit')
parser.add_argument('--expected_sskingdom', help='Expected sskingdom value for this dataset')
parser.add_argument('--evalue', help='evalue for blastn search')
args = parser.parse_args()

# search sequences in nt database
ret = call("blastn -db nt -outfmt '6 delim=@#@# 'qaccver saccver sskingdom ssciname salltitles pident evalue' -num_threads " + args.num_threads +
	" -out blastn_report.csv -query " + args.fasta + "-evalue " + args.evalue, shell=True)
if ret!=0:
    print("blastn search failed")
    exit()

nt = pd.read_csv("blastn_report.csv", sep='@#@#',header=None, engine='python')
nt.columns = ['qaccver', 'saccver', 'sskingdom', 'ssciname', 'salltitles', 'pident', 'evalue']

# select best hits
nt = nt.sort_values(by=['evalue'], ascending=True).drop_duplicates(['qaccver'])
ids_to_remove=nt[nt['sskingdom']!=args.expected_sskingdom]['qaccver']

ids_to_remove.to_csv("toremovelist", header = False, index = False)

# remove sequences that were assigned to the wrong kingdom or don'have sskingdom value
if not os.path.isfile("toremovelist") or os.path.getsize("toremovelist") == 0:
	ret = call('cp '+ args.fasta  + ' filtered.fasta', shell=True)
	if ret!=0:
		print("Copy file failed")
		exit()
else:
	ret = call('seqkit grep -j '+args.num_threads+' -v -f toremovelist ' + args.fasta + ' > filtered.fasta', shell=True)
	if ret!=0:
	    print("Seqkit grep failed")
	    exit()