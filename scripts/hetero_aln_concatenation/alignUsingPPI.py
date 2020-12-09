#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 01:01:58 2020

@author: farhan
"""

import os, sys

def readPPIDict(file):
    ppi_dict={}
    with open (file) as f:
        for line in f:
            line=line.strip()
            key=line.split(":")[0].strip()
            values=line.split(":")[1].strip().split(",")
            ppi_dict[key]=values
    return ppi_dict

def readIndexFile(idx_file):
    idx_list=[]
    with open (idx_file) as f:
        for line in f:
            idx_list.append(line.strip())
    return idx_list

ppi_dict_file=os.path.abspath(sys.argv[1])
idx_file_A=os.path.abspath(sys.argv[2])
idx_file_B=os.path.abspath(sys.argv[3])

aln_list=[]
ppi_dict=readPPIDict(ppi_dict_file)
P=readIndexFile(idx_file_A)
Q=readIndexFile(idx_file_B)
lenP=len(P)
lenQ=len(Q)
