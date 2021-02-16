#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 15:07:20 2020

@author: farhan
"""

import os, sys
from DNCON_lib import *
#from readRR import readRRFile
import numpy as np

max_len=400
path_of_lists="/data/farhan/SoftwareTools/docking_outputs/precision_code_heteromer/"
test_dict_L=build_dataset_dictionaries_file(path_of_lists,"raj_test_list.txt")
print (len(test_dict_L))
#valid_dict_L,_=build_dataset_dictionaries_test(path_of_lists)
#train_dict_L,_=build_dataset_dictionaries_train(path_of_lists)
test_dict=subset_pdb_dict(test_dict_L,30,max_len,10000,"ordered")
#valid_dict=subset_pdb_dict(valid_dict_L,30,max_len,10000,"ordered")
#train_dict=subset_pdb_dict(train_dict_L,30,max_len,10000,"ordered")


test_name_dict=list(test_dict.keys())
#valid_name_dict=list(valid_dict.keys())
#train_name_dict=list(train_dict.keys())

rr_dir="/data/farhan/SoftwareTools/docking_outputs/test_raj_400/"#"/data/farhan/SoftwareTools/DeepComplex/training/dncon4_net_hetero/models/custom/DNCON4_RESOTHER_HETER30_otherfea_filter64_layers34_ftsize3_binary_crossentropy_3.0/predict_map/test/"
label_dir="/data/farhan/SoftwareTools/docking_outputs/test_raj_400/"#"/data/farhan/SoftwareTools/DeepComplex/training/dncon4_net_hetero/models/custom/DNCON4_RESOTHER_HETER30_otherfea_filter64_layers34_ftsize3_binary_crossentropy_3.0/predict_map/test/" #"/data/farhan/SoftwareTools/HomopolymerProject/data/homodimers/scripts/Y-Labels/"
valid_outdir="/data/farhan/SoftwareTools/DeepComplex/training/dncon4_net_hetero/models/custom/DNCON4_RESOTHER_HETER30_otherfea_filter64_layers34_ftsize3_binary_crossentropy_3.0/predict_map/val_prec/"
test_outdir="/data/farhan/SoftwareTools/DeepComplex/training/dncon4_net_hetero/models/custom/DNCON4_RESOTHER_HETER30_otherfea_filter64_layers34_ftsize3_binary_crossentropy_3.0/predict_map/test_prec_raj_400/"#"/data/farhan/SoftwareTools/DeepComplex/training/dncon4_net_intra/models/custom/DNCON4_RESOTHER_DEEPCOV_otherfea_filter64_layers34_ftsize3_4.0/predict_map/precision_test_zdock/"
train_outdir="/data/farhan/SoftwareTools/DeepComplex/training/dncon4_net_hetero/models/custom/DNCON4_RESOTHER_HETER30_otherfea_filter64_layers34_ftsize3_binary_crossentropy_3.0/predict_map/train_prec/"
#key="2NLF"
#os.system("python getPrecision_inter_twocmaps_v2.py "+rr_dir+"Y-"+key+".txt "+rr_dir+key+"_.rr > "+valid_outdir+key+"_inter_precision.txt")
i=0
print (len(test_name_dict))
#sys.exit()
"""
for key in train_name_dict:
    os.system("python getPrecision_inter_twocmaps_v2.py "+rr_dir+"Y-"+key+".txt "+rr_dir+key+"_.rr > "+train_outdir+key+"_inter_precision.txt")
"""

"""
for key in valid_name_dict:
    os.system("python getPrecision_inter_twocmaps_v2.py "+rr_dir+"Y-"+key+".txt "+rr_dir+key+"_.rr > "+valid_outdir+key+"_inter_precision.txt")
    #break

"""

for key in test_name_dict:
    i+=1
    if not (os.path.isdir(test_outdir)): os.makedirs(test_outdir)
    if not (os.path.exists(rr_dir+key+"_.rr")): continue
    print (str(i)+" Calculating precision for: "+key)
    os.system("python getPrecision_inter_twocmaps_hetero.py "+label_dir+"Y-"+key+".txt "+rr_dir+key+"_.rr > "+test_outdir+key+"_inter_precision.txt")
    #break

"""
with open ("test.txt","w") as f:
    for key in test_name_dict:
        f.write(key+"\n")
with open ("valid.txt","w") as f:
    for key in valid_name_dict:
        f.write(key+"\n")

with open ("train.txt","w") as f:
    for key in train_name_dict:
        f.write(key+"\n")
"""

