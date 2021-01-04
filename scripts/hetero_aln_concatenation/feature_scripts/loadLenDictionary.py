#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 16:00:58 2020

@author: farhan
"""

#this script loads the fasta_length_dictionary.txt file into a dictionary
#usage: python loadLenDictionary.txt <fasta_dictionary.txt>

def loadLenDictionary(dict_file):
    fasta_dict={}
#    i=0
    with open(dict_file,"r") as f:
        for line in f:
#            i+=1
            fasta_dict[line.strip().split()[0].strip()]=int(line.strip().split()[1].strip())

#    print (len(fasta_dict.keys()),i)
    return fasta_dict
