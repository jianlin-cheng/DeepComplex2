print('here')
import os,sys

from readPDBColumns_hetero import readPDB,readAtom,write2File,contents2Info,addColumn,reassembleLines,getChain,getName
import numpy as np
#python pdb2distance_inter_heavy_hetero.py /home/rajroy/test/real_pdb/1AVAA.pdb /home/rajroy/test/real_pdb/1AVAC.pdb 6.0 /home/rajroy/test/1AVA

def pdb2FastaFromSplitContents(split_contents):
    fst=""
    letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V'}
    prev_res_num=split_contents[0]["res_num"]
    fst+=letters[split_contents[0]["res_name"]]
    for items in split_contents:
        if (items["res_num"]==prev_res_num): continue
        fst+=letters[items["res_name"]]
        prev_res_num=items["res_num"]

    return fst

def getCoordinate(atom):
    coordinate={}
    coordinate["x"]=atom["x"]
    coordinate["y"]=atom["y"]
    coordinate["z"]=atom["z"]
    return coordinate

def distance(coord1, coord2):
    x1=float(coord1["x"])
    y1=float(coord1["y"])
    z1=float(coord1["z"])
    x2=float(coord2["x"])
    y2=float(coord2["y"])
    z2=float(coord2["z"])
    d=np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

    return d

def createDistanceMapAllAtoms(atom_list_A,atom_list_B,dist):
    result_list_string=[]

    for atom_A in atom_list_A:
        for atom_B in atom_list_B:
            string=""
            if (atom_A["chain"]==" " or atom_A["chain"]=="  " or atom_A["chain"]==""): atom_A["chain"]=chain_1
            if (atom_B["chain"]==" " or atom_B["chain"]=="  " or atom_B["chain"]==""): atom_B["chain"]=chain_2

            string+=atom_A["chain"]+" "+atom_A["serial"]+" "+atom_A["res_name"]+" "+atom_A["res_num"]+" "+atom_A["x"]+" "+atom_A["y"]+" "+atom_A["z"]+" "+atom_A["atom_name"]+" | "+ atom_B["atom_name"]+" "+atom_B["chain"]+" "+atom_B["serial"]+" "+atom_B["res_name"]+" "+atom_B["res_num"]+" "+atom_B["x"]+" "+atom_B["y"]+" "+atom_B["z"]+" | " + str(distance(getCoordinate(atom_A),getCoordinate(atom_B)))


            if (distance(getCoordinate(atom_A),getCoordinate(atom_B))<dist): result_list_string.append(string+"\n")

    return result_list_string

def removeRedundantContacts(contact_list):
    new_list=[]
    contact_dict={}

    for contact in contact_list:
        x_y=contact.split()[0]+" "+contact.split()[1]
        if x_y in contact_dict:
            temp_array=contact_dict[x_y]
        else:
            temp_array=[]
        temp_array.append(float(contact.split()[4]))
        contact_dict[x_y]=temp_array
    # i j 0 d 1/d
    for key in contact_dict.keys():
        min_value= min(contact_dict[key])
        temp_string =  key+' 0 '+str(min_value)[0:5]+ ' '+str(1/ min_value)[0:5]+"\n"
        new_list.append(temp_string)

    return new_list


def createDistanceMapHeavyAtoms(atom_list_A,atom_list_B,dist=6.0): #new version. Chose all atoms but the hydrogens
    result_list_string=[]
    distance_list_string=[]
    rr_list_string=[]
    for atom_A in atom_list_A:
        for atom_B in atom_list_B:
            string=""
            if (atom_A["chain"]==" " or atom_A["chain"]=="  " or atom_A["chain"]==""): atom_A["chain"]=chain_1
            if (atom_B["chain"]==" " or atom_B["chain"]=="  " or atom_B["chain"]==""): atom_B["chain"]=chain_2

            temp_dist= 0.0
            string=atom_A["chain"]+" "+atom_A["serial"]+" "+atom_A["res_name"]+" "+atom_A["res_num"]+" "+atom_A["x"]+" "+atom_A["y"]+" "+atom_A["z"]+" "+atom_A["atom_name"]+" | "+ atom_B["atom_name"]+" "+atom_B["chain"]+" "+atom_B["serial"]+" "+atom_B["res_name"]+" "+atom_B["res_num"]+" "+atom_B["x"]+" "+atom_B["y"]+" "+atom_B["z"]+" | " + str(distance(getCoordinate(atom_A),getCoordinate(atom_B)))
            if (atom_A["element"].strip()=="H" or atom_A["element"].strip()=="D"): continue
            if (atom_B["element"].strip()=="H" or atom_A["element"].strip()=="D"): continue
            temp_dist=distance(getCoordinate(atom_A),getCoordinate(atom_B))
            # if (temp_dist<=dist): 
            # if (distance(getCoordinate(atom_A),getCoordinate(atom_B))<=dist): 
            distance_string=atom_A["res_num"].strip()+" "+atom_B["res_num"].strip()+" 0 "+str(dist)+" "+ str(temp_dist)[0:5]
            # rr_string=atom_A["res_num"].strip()+" "+atom_B["res_num"].strip()+" 0 "+str(temp_dist)[0:5]+" "+ str(0.5+1/(temp_dist))
            result_list_string.append(string+"\n")
            distance_list_string.append(distance_string+"\n")
            # rr_list_string.append(rr_string+"\n")

    distance_list_string=removeRedundantContacts(distance_list_string)
    # rr_list_string=removeRedundantContacts(rr_list_string)
    return result_list_string,distance_list_string,rr_list_string

def writeToFile(outfile,stuff,fasta_1,fasta_2,name_A,name_B):
    with open (outfile,"w") as f:
        f.write("#"+name_A+name_B+"\n")
        f.write(fasta_1+"\n")
        f.write(fasta_2+"\n")
        f.writelines(stuff)
    return

def readFastaDict(fasta_dict_file="fasta_dictionary.txt"):
    fasta_dict={}
    if not (os.path.exists(fasta_dict_file)):
        print(fasta_dict_file+" fasta dictionary file not found. Exiting")
        sys.exit(-6)
    with open (fasta_dict_file,"r") as f:
        for line in f:
            name=line.strip().split(":")[0].strip()
            fast=line.strip().split(":")[1].strip()
            fasta_dict[name]=fast
    return fasta_dict

def getFastaFromDictionary(pdbfile,fasta_dict):
    #fasta_dict=readFastaDict(fasta_dict_file)
    name=pdbfile.split("/")[-1].replace(".pdb","").replace(".fasta","").replace(".atom","")
    return fasta_dict[name]

# TEST_CASE
# pdbfile_A="/home/rajroy/test/real_pdb/1AVAA.pdb"
# pdbfile_B="/home/rajroy/test/real_pdb/1AVAC.pdb"
# dist="6.0"
# outfile="/home/rajroy/test/1AVA_out"


pdbfile_A=sys.argv[1]
pdbfile_B=sys.argv[2]
dist=sys.argv[3]
outfile=sys.argv[4]
pdb=getName(pdbfile_A)

#set the fasta dictionary here:
fasta_dict_file="/home/rajroy/fasta_dictionary.txt"
fasta_dict=readFastaDict(fasta_dict_file)

split_contents_A=contents2Info(readPDB(pdbfile_A))
split_contents_B=contents2Info(readPDB(pdbfile_B))
chain_1=getChain(pdbfile_A)
chain_2=getChain(pdbfile_B)
name_A=pdb+chain_1
name_B=pdb+chain_2
atom_list_A=split_contents_A
atom_list_B=split_contents_B
#Get fasta from fasta_dictionary.txt
#fasta_A=pdb2FastaFromSplitContents(atom_list_A)
fasta_A=getFastaFromDictionary(pdbfile_A,fasta_dict)
L_A=len(fasta_A)
#fasta_B=pdb2FastaFromSplitContents(atom_list_B)
fasta_B=getFastaFromDictionary(pdbfile_B,fasta_dict)
L_B=len(fasta_B)

if (fasta_A=="" or fasta_B==""):
    print("No fasta for "+pdb)
    os.system("echo '"+pdb+"' >> no_fasta.txt")
    sys.exit("No fasta found for "+pdb)

result_list,pw_dist,pw_rr=createDistanceMapHeavyAtoms(atom_list_A,atom_list_B,float(dist))
print (len(pw_dist),len(pw_rr), len(result_list))
if (len(pw_dist)==0):
    print (pdb+" "+": No interchain contacts less than "+dist+"...")
    os.system("echo "+pdb+" "+": No interchain contacts less than "+dist+" ... >> not_done_reason_heavy.txt")
    sys.exit("-2")

outfile_dist=outfile.replace(".txt","")+"_dist_"+".txt"
outfile=outfile+"_"+chain_1+chain_2+".txt"

write_flag=True

if (write_flag):
    writeToFile(outfile_dist,pw_dist,fasta_A,fasta_B,name_A,name_B)
    # writeToFile(outfile.replace(".txt",".rr"),pw_rr,fasta_A,fasta_B,name_A,name_B)


if not(os.path.exists(outfile_dist)):
    sys.exit("-4")



