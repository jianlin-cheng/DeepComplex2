#this scripts reads the binding file and creates a PPI dictionary of indexes.
import sys

binding_file = sys.argv[1]
binding_dict={}
binding_pair_list=[]
binding_starter_list=[]

with open (binding_file) as f:
    for line in f:
        p=line.strip().split()[0].strip()
        q=line.strip().split()[1].strip()
        binding_pair_list.append([p,q])
        binding_starter_list.append(p)

for x in binding_starter_list:
    cslist=""
    for pair in binding_pair_list:
        if x==pair[0]:
            #print (pair)
            cslist+=pair[1]+","
    cslist=cslist.rstrip(",")
    binding_dict[x]=cslist

with open ("ppi_dict.txt","w") as f:
    for key, value in binding_dict.items():
        f.write(key+":"+value+"\n")
        
