#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 02:46:52 2021

@author: farhan
"""
import os, sys
import numpy as np
from loadFastaDictionary import loadFastaDictionary
from loadLenDictionary import loadLenDictionary

def getNewLabelFromOld(filename,len_A, len_B):
    old_Y=np.loadtxt(filename)
    print("len_A=",len_A,"len_B=",len_B)
    all_ones=np.sum(old_Y)
    slice_Y=old_Y[0:len_A,len_A:]
    slice_ones=np.sum(slice_Y)
    print ("all_ones=",all_ones,": slice_ones=",slice_ones)
    print("Old_shape=",old_Y.shape)
    print("Slice_Y_shape=",slice_Y.shape)
    if (all_ones/2 != slice_ones): 
        mismatch_list.append(filename+"\n")
    return slice_Y

mismatch_list=[]
old_label_dir="/data/farhan/SoftwareTools/DeepComplex/training/dncon4_net_hetero/data/HETER30/Y_label/"
new_label_dir="/data/farhan/SoftwareTools/DeepComplex/training/dncon4_net_hetero/data/HETER30/Y-Label-reduced/"

fasta_dict=loadFastaDictionary("all_hetero_monomer_fasta_dict.txt")
len_dict=loadLenDictionary("all_hetero_monomer_fasta_len.txt")
all_heteromer_list=[]
with open ("all_heterodimer_list.txt") as f:
    for line in f:
        all_heteromer_list.append(line.strip())

for dimer in all_heteromer_list:
    print (dimer)
    if "__" in dimer:
        monomer_A=dimer.split("__")[1]
        monomer_B=dimer.split("__")[0]
        Y_file=monomer_A+"_"+monomer_B+".txt"
        new_Y_file=monomer_B+"_"+monomer_A+".txt"
        new_label=getNewLabelFromOld(old_label_dir+Y_file,len_dict[monomer_A],len_dict[monomer_B])
        #print (new_label.shape)
        #print (len_dict[monomer_A],len_dict[monomer_B])
        slice_Y=np.transpose(new_label)
        np.savetxt(new_label_dir+new_Y_file,slice_Y,fmt="%.1f")
        np.savetxt(new_label_dir+Y_file,new_label,fmt="%.1f")
        #print ("new_shape=",slice_Y.shape)
        #print (len_dict[monomer_B],len_dict[monomer_A])
        continue
    if "_" in dimer:
        monomer_A=dimer.split("_")[0]
        monomer_B=dimer.split("_")[1]
        Y_file=monomer_A+"_"+monomer_B+".txt"
        new_Y_file=monomer_A+"_"+monomer_B+".txt"
        new_label=getNewLabelFromOld(old_label_dir+Y_file,len_dict[monomer_A],len_dict[monomer_B])
        np.savetxt(new_label_dir+Y_file,new_label,fmt="%.1f")
        #print ("extracted_shape=",new_label.shape)
        #print (len_dict[monomer_A],len_dict[monomer_B])
        #slice_Y=np.transpose(new_label)
        #print ("new_shape=",slice_Y.shape)
        #print (len_dict[monomer_B],len_dict[monomer_A])
        #continue
    #break

with open ("mismatch_list.txt","w") as f:
    f.writelines(mismatch_list)

