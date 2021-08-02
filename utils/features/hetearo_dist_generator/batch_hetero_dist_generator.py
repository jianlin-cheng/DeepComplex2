import os


def file_reader(_input_dir):
    f = open(_input_dir, "r")
    if f.mode == 'r':
        contents = f.read().splitlines()
        f.close()
    out_array=[]
    for val in contents:
        final_name = ""
        if "__" in val:
            temp = val.replace("__","_")
            final_name= temp.split("_")[1]+'_' +temp.split("_")[0]
        else:
            final_name = val.split("_")[0] + '_' + val.split("_")[1]
        out_array.append(final_name)
    return out_array


def loadFastaDictionary(dict_file):
    fasta_dict = {}
    with open(dict_file, "r") as f:
        for line in f:
            fasta_dict[line.strip().split(":")[0].strip()] = line.strip().split(":")[1].strip()
    return fasta_dict

fasta_dict = loadFastaDictionary('/home/rajroy/fasta_dictionary.txt')

prog_dir = "/home/rajroy/Documents/DeepComplex2/utils/features/hetearo_dist_generator/pdb2distance_inter_heavy_hetero.py"
pdb_file_dir = "/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/rajroy/pdb_atom/"
input_list = "/home/rajroy/het30_dup_list.txt"
output_dir = "/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/rajroy/NEW_DIST_FILES/Dist_files/"
input_list_array=file_reader(input_list)
input_list_array=list(dict.fromkeys(input_list_array))

for pdb in input_list_array:

    name_array=pdb.split("_")
    total =len(fasta_dict.get(name_array[0]))+len(fasta_dict.get(name_array[1]))
    if total <601:
        cmd = "python "+prog_dir +" "+pdb_file_dir+name_array[0]+".atom"+" "+pdb_file_dir+name_array[1]+".atom" +" 6.0 "+ output_dir+pdb
        print(cmd)
        os.system(cmd)
