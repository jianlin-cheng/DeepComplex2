def search(d,v):
    for k in d:
        if k.strip()==v.strip():
            return True
    return False

dict_keys=[]
aln_keys=[]

with open ("ppi_dict_keys.txt") as f:
    for line in f:
        dict_keys.append(line.strip())

with open ("missing.txt") as f:
    for line in f:
        aln_keys.append(line.strip())


print (len(aln_keys))
print (len(dict_keys))

aln_set=set(aln_keys)
dict_set=set(dict_keys)

print (len(dict_set-aln_set))
value=""
search(dict_key,value)