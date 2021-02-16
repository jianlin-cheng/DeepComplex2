#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:01:04 2019

@author: farhan
"""
#modified for heteromers
#sorts an rr file and prints its contents
#usage: python sortrr.py <contact_map_file.rr> <reverse=True (high prob to low) or False> 
import os,sys

def prob(e):
    return e[-1]


rrfile=sys.argv[1]
rev=sys.argv[2]
if not os.path.exists(rrfile): sys.exit(rrfile+" not found. Quitting!")
def readFile(file):
    contents=[]
    fasta_A=""
    fasta_B=""
    line=""
    lab=""
    with open (file,"r") as f:
        ln=f.readline().strip()
        if (ln.startswith("#")): 
            lab=ln.strip()
        fasta_A=f.readline().strip()
        fasta_B=f.readline().strip()
        #print (lab)
        #print (fasta_A)
        #print (fasta_B)
        for line in f:
            if (line.strip().startswith("0") or line.strip().startswith("1") or line.strip().startswith("2") or line.strip().startswith("3") or line.strip().startswith("4") or line.strip().startswith("5") or line.strip().startswith("6") or line.strip().startswith("7") or line.strip().startswith("8") or line.strip().startswith("9")):
                contents.append(line.strip())
                continue
            if (line.strip().startswith("#") or line.strip().startswith("PFRMAT") or line.strip().startswith("METHOD") or line.strip().startswith("MODEL") or line.strip().startswith("REMARK") or line.strip().startswith("TARGET") or line.strip().startswith("AUTHOR")):
                continue
#            fasta_A=line.strip()
#    if (line.strip()=="END"): fasta=fasta[0:len(fasta)-3]
    return lab,fasta_A,fasta_B,contents

label,fasta_A,fasta_B,contents=readFile(rrfile)

def toNum(contents):
    new_contents=[]
    for line in contents:
        split=line.split()
        i=int(split[0])
        j=int(split[1])
        zero=int(split[2])
        dist=int(float(split[3]))
        val=float(split[4])
        #print (i,j,zero,dist,val)
        new_contents.append([i,j,zero,dist,val])
    return new_contents



contents=toNum(contents)
if rev.strip()=="True": contents.sort(key=prob,reverse=True)
if rev.strip()=="False": contents.sort(key=prob,reverse=False)

#print (contents)
print (label)
print(fasta_A)
print(fasta_B)
for item in contents:
    print(str(item[0])+" "+str(item[1])+" "+str(item[2])+" "+str(item[3])+" "+str(item[4]))#+" "+str(item[5]))
#print("END")

