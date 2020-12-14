s=""
with open ("ppi_dict_temp.txt") as f:
    for line in f:
        print (line.strip()[0:100])