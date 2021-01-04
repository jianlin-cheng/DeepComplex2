#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 02:03:39 2021

@author: farhan
"""

import os, sys
from glob import glob

fasta_dir=os.path.abspath(sys.argv[1])+"/"
all_files=glob(fasta_dir+"*.fasta")

print (len(all_files))
#print (os.path.basename(all_files[0]))
#print (all_files[1])
#print (all_files[2])
#print (all_files[3])
#print (os.path.basename(all_files[4]))
#sys.exit()

len_dict={}

for file in all_files:
    file_name=os.path.basename(file)
    if "_" in file_name or "__" in file_name: continue
    fasta=""    
    with open (file) as f:
        for line in f:
            if line.strip().startswith(">"): continue
            fasta+=line.strip()
    len_dict[file_name]=str(len(fasta))
with open ("all_hetero_monomer_fasta_len.txt","w") as f:
    for key, value in len_dict.items():
        f.write(key+"\t"+value+"\n")

