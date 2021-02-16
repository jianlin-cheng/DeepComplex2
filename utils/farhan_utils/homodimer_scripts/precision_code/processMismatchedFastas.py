#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 17:32:42 2020

@author: farhan
"""

#this script processes the mismatched fasta sequence during evaluation of two homodimer contact map (.rr) files:
#1. Searching the fasta_dictionary.txt to find the fasta sequences of the respective .atom files
#2. Aligns the sequences and finds sequence similarity using findSequenceSimilarity.py script
#3. Reindexes the respective .atom files to a selected destination "reindexed_atom"
#Usage: python processMismatchedFastas.py <PDB_NAME> <cmap_file1> <cmap_file2> <outfolder>
import os,sys
from readRR import readRRFile

def readFastaDict(fasta_dict_file):
    fasta_dict={}
    with open (fasta_dict_file,"r") as f:
        for line in f:
            fasta_dict[line.strip().split(":")[0].strip()]=line.strip().split(":")[1].strip()
    return fasta_dict


pdb_name=sys.argv[1]#"2HXR"#"2EHP"#"4I24"#"1BFT"#"1UZY"#"5IZV"#"11AS"
rr_file1=os.path.abspath(sys.argv[2])
rr_file2=os.path.abspath(sys.argv[3])
outdir=os.path.abspath(sys.argv[4])+"/"
if not os.path.isdir(outdir):os.makedirs(outdir)
print("Processing fasta mismatches in chains for: "+pdb_name)
print("Reading fasta_dictionary...")

fasta1,rr1=readRRFile(rr_file1)
fasta2,rr2=readRRFile(rr_file2)
name1=os.path.basename(rr_file1).split(".")[0]
name2=os.path.basename(rr_file2).split(".")[0]
fasta_dict={name1:fasta1,name2:fasta2}#readFastaDict("fasta_dictionary.txt") #reads the fasta_dictionary.txt file as a dictionary
print (fasta_dict)
with open ("rr_fasta_dict.txt","w") as f:
    for key, value in fasta_dict.items():
        f.write(key+":"+value+"\n")
key_list=list(fasta_dict.keys()) #get the keys of the dictionary
this_key_list=[] #stores the list of keys for pdb_name
print("Done!")

for key in key_list:
    if pdb_name in key:
        this_key_list.append(key)
print ("This key_list:",this_key_list)
print ("Atom file list found: ",this_key_list)
print("key_list:",fasta_dict.keys())
print(fasta_dict[this_key_list[0]])
print(fasta_dict[this_key_list[1]])
#the following arranges the fasta lengthwise. x=lower length, y= higher length
if len(fasta_dict[this_key_list[0]]) > len(fasta_dict[this_key_list[1]]):
    y=fasta_dict[this_key_list[0]]
    key_y=this_key_list[0]
    x=fasta_dict[this_key_list[1]]
    key_x=this_key_list[1]
    print("y > x", len(y),len(x))
else:
    y=fasta_dict[this_key_list[1]]
    key_y=this_key_list[1]
    x=fasta_dict[this_key_list[0]]
    key_x=this_key_list[0]
    print("x > y", len(x),len(y))
#No need if chain sequences are similary
if x==y: sys.exit("All fasta sequences are similar. Quitting!")
#2. and #3.
exit_code=os.system("python findSequenceSimilarity.py "+x+" "+y+" "+pdb_name+" "+key_x+" "+key_y)
if (exit_code!=0):
    sys.exit("Some thing went wrong for "+pdb_name+". Quitting!")

#Now reindex the .rr files
map_file1="./reindexed_mapping_function/"+key_x+"_map.txt"
seq_aln_file="./aligned_seq_folder/"+pdb_name+".aln.txt"
outfolder=outdir
print(key_x,rr_file2,map_file1)
os.system("python reindex_rr.py "+rr_file2+" "+map_file1+" "+seq_aln_file+" "+outfolder+" "+key_x)

#Fix the sequence problem later

"""
if key_x in rr_file1 and key_x in map_file1:
    os.system("python reindex_rr.py "+rr_file1+" "+map_file1+" "+seq_aln_file+" "+outfolder+" "+key_x)
else:
    os.system("python reindex_rr.py "+rr_file2+" "+map_file1+" "+seq_aln_file+" "+outfolder+" "+key_x)
"""

map_file2="./reindexed_mapping_function/"+key_y+"_map.txt"
seq_aln_file="./aligned_seq_folder/"+pdb_name+".aln.txt"
outfolder=outdir
print(key_y,rr_file1,map_file2)
#os.system("python reindex_rr.py "+rr_file2+" "+map_file2+" "+seq_aln_file+" "+outfolder+" "+key_y)
os.system("python reindex_rr.py "+rr_file1+" "+map_file2+" "+seq_aln_file+" "+outfolder+" "+key_y)

