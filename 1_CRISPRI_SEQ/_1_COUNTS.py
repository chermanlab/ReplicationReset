import os
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import gzip
from tqdm import tqdm
import pandas as pd

#We first open the library file as it was synthesized to obtain the sequence of all guides.
library_handle=open('_library.txt','r')
library=SeqIO.parse(library_handle,'fasta')

os.chdir('./0_FORSRA')
files=os.listdir(os.getcwd())

#We then create a dictionary to store count data: 
#It contains 1 key per guide, then each value is a dictionary with 1 key and value per sample.
targets={}
for l in library:
    targets[str(l.seq[:20])]={}
    for file in files:
        sample=file.split('.fastq')[0]
        targets[str(l.seq[:20])][sample]=0

#We then go through each fastq sequencing file (1 per condition and replicate).
#If the read is in the library, then we store it in the dictionary.
output_file = "./1.1_count_stats.txt"
sys.stdout = open(output_file, "w")

for file in tqdm(files):
    correct_sequences=0
    total_reads=0
    sample=file.split('.fastq')[0]
    with open(file,'rt') as handle:
        for r in SeqIO.parse(handle,'fastq'):
            total_reads+=1
            if str(r.seq[:20]) in targets.keys():
                targets[r.seq[:20]][sample]+=1
                correct_sequences+=1
    print(sample+":")
    print(str(total_reads)+' reads')
    print(str(round(100*correct_sequences/total_reads,2))+"% correct reads")

sys.stdout.close()
sys.stdout = sys.__stdout__

#The counts table in finally saved as a dataframe.
counts_table=pd.DataFrame.from_dict(targets, orient='index')
os.chdir('./')
counts_table.to_csv('./1.0_count_table.txt',sep='\t')
