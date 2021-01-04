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

#print (len(all_files))
#print (os.path.basename(all_files[0]))
#print (all_files[1])
#print (all_files[2])
#print (all_files[3])
#print (os.path.basename(all_files[4]))
#sys.exit()

dimer_list=[]

for file in all_files:
    file_name=os.path.basename(file)
    if not ("_" in file_name or "__" in file_name): 
        #all_files.remove(file)
        continue
    file_name=file_name.replace(".fasta","")
    if "__" in file_name:
        name_A=file_name.split("__")[0]
        name_B=file_name.split("__")[1].split("_")[0]
        file_name=name_A+"__"+name_B
        dimer_list.append(file_name+"\n")
        continue
    if "_" in file_name:
        name_A=file_name.split("_")[0]
        name_B=file_name.split("_")[1]
        file_name=name_A+"_"+name_B
    
    dimer_list.append(file_name+"\n")
    
    
with open ("all_heterodimer_list.txt","w") as f:
    f.writelines(dimer_list)

