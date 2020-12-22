import os, sys
from glob import glob

alnfolder=os.path.abspath(sys.argv[1])
outfile=sys.argv[2]

if os.path.exists(outfile):os.remove(outfile)
alnlist=glob(alnfolder+"/*.aln")
#print (len(alnlist))
#print (alnfolder)
for file in alnlist:
#    print (file)
#    break
    os.system("python calc_Neff.py "+file.strip()+" "+outfile)

"""
with open("aln_list.lst","r") as f:
    for alnfile in f:
        os.system("echo "+alnfile.strip().replace(".aln","")+" >> tmp.eff && python calc_Neff.py "+alnfile.strip()+" >> tmp.eff")

l=[]
with open ("tmp.eff","r") as f:
    for line in f:
        s=""
        if (line.strip().startswith("T")):
            s+=line.strip()
            line=f.readline()
            s+="\t"+line.strip()
            l.append(s+"\n")

l[len(l)-1]=l[len(l)-1].strip()

os.remove("tmp.eff")

with open("N_eff.txt","w")as f:
    for items in l:
        f.write(items)
"""
