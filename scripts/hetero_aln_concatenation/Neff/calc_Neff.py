import os,sys,shutil

def getName(file):
    if ("/" in file):
        split = file.split("/")
        file = split[len(split)-1]
    
    return file.replace(".aln","")

alnstat="/storage/htc/bdm/tools/MetaPSICOV/bin/alnstats"
alnfile=sys.argv[1]
outfile=sys.argv[2]
L=0
line_count=0
with open (alnfile) as f:
    for line in f:
        L=len(line.strip())
#        break
        line_count+=1

if (not os.path.isdir(os.getcwd()+"/temp_alnstat")):
    os.mkdir(os.getcwd()+"/temp_alnstat")

pairstatsfile=os.getcwd()+"/temp_alnstat/"+getName(alnfile)+".pairstats"
colstatsfile=os.getcwd()+"/temp_alnstat/"+getName(alnfile)+".colstats"
#print(alnstat+" "+alnfile+" "+colstatsfile+" "+pairstatsfile)
os.system(alnstat+" "+alnfile+" "+colstatsfile+" "+pairstatsfile)
with open(colstatsfile,"r") as f:
    for line in f:
        line=f.readline()
        Neff=f.readline().strip()
        break
print(Neff)
with open (outfile,"a+") as f:
    f.write(getName(alnfile)+","+str(L)+","+str(line_count)+","+str(Neff)+"\n")
shutil.rmtree(os.getcwd()+"/temp_alnstat")    
