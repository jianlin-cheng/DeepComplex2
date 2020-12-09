#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 01:02:52 2020

@author: farhan
"""
#this script will create a fasta dictionary and and index file using the a3m output file
import os,sys

a3m_file=os.path.abspath(sys.argv[1])
dirname=os.path.dirname(a3m_file)+"/"
name=os.path.basename(a3m_file).split("_")[0].replace(".fasta","")
index_list=[]
aln_dict={}

with open (a3m_file,"r") as f:
    for line in f:
        if line.startswith(">"):
            idx=line.strip().replace(">","").split()[0] #take just the id of the index
            index_list.append(idx+"\n")
            line=f.readline().strip()
            aln_dict(index_list[-1].strip())=line.strip()
index_list[-1]=index_list[-1].strip()

with open (dirname+name+".idx","w") as f:
    f.writelines(index_list)

with open (dirname+name+"_aln_dict.txt","w") as f:
    for key, value in aln_dict.items():
        f.write(key+":"+value+"\n")
        