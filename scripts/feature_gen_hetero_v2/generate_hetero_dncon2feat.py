#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 02:21:47 2020

@author: farhan
"""
import os, sys
from glob import glob

#this script will generate the dncon2 features from a list of files. The list needs to be is format: FASTAA_FASTAB
#Also supply the alignmet folder and the fasta folder
usage=sys.argv[0]+" <list_file> <fasta_folder> <alignment_folder> <outfolder>"

list_file=os.path.abspath(sys.argv[1])
fasta_folder=os.path.abspath(sys.argv[2])+"/"
aln_folder=os.path.abspath(sys.argv[3])+"/"
outdir=os.path.abspath(sys.argv[4])+"/"
paths_file="/storage/htc/bdm/farhan/string_db/hetero30_aln_run/feature_gen_hetero_v2/paths_lewis.txt"
if len(sys.argv)!=5: 
    print("Needs 4 input parameters")
    sys.exit(usage)

if not os.path.exists(list_file): sys.exit(list_file+" not found. Quitting!")
if not os.path.isdir(fasta_folder): sys.exit(fasta_folder+" not found. Quitting!")
if not os.path.isdir(aln_folder): sys.exit(aln_folder+" not found. Quitting!")
if not os.path.isdir(outdir): os.makedirs(outdir)

file_list=[]
with open (list_file) as f:
    for line in f:
        if "__" in line: line=line.replace("__","_")
        line=line.strip().split("_")[0]+"_"+line.strip().split("_")[1]
        name_A=line.split("_")[0]
        name_B=line.split("_")[1]
        if not os.path.exists(fasta_folder+name_A+".fasta"): continue
        if not os.path.exists(fasta_folder+name_B+".fasta"): continue
        if not os.path.exists(aln_folder+name_A+"_"+name_B+".a3m"): continue
        if not os.path.exists(aln_folder+name_A+"_"+name_B+".aln"): 
            os.system("grep -v '^>'"+ aln_folder+name_A+"_"+name_B+".a3m | sed 's/[a-z]//g' > "+aln_folder+name_A+"_"+name_B+".aln")
        file_list.append(line.strip())

for file in file_list:
    fasta_A=file.split("_")[0]+".fasta"
    fasta_B=file.split("_")[1]+".fasta"
    os.system("sbatch /storage/htc/bdm/farhan/string_db/hetero30_aln_run/dncon2_feat_run/generate_dncon2_feat_hetero.sh "+fasta_folder+fasta_A+" "+
              fasta_folder+fasta_B+" "+outdir+file+" "+paths_file)


