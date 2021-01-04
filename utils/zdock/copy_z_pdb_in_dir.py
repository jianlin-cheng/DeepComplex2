import os


def single_fasta_file_reader(_seq_file):
    file = open(_seq_file, "r")
    output_array = []
    if file.mode == 'r':
        output_array = file.read().splitlines()
        file.close()

    return output_array


input_list = "/home/raj/test_list.txt"
fasta_files = single_fasta_file_reader(input_list)
#ginger
pdb_dir = '/home/casp14/Data/pdb_no_blanks_alt_removed/'
out_dir = "/home/raj/pdb_test_400/"
missing_count= 0
for pdb in fasta_files:
    name = ""
    if '__' in pdb:
        name = pdb.replace('__', '_')
    else:
        name =pdb
    name_array = name.split('_')
    if not os.path.exists(pdb_dir+name_array[0]+'.atom'):
        missing_count = missing_count+1
    if not os.path.exists(pdb_dir+name_array[1]+'.atom'):
        missing_count = missing_count + 1
    cmd_1 = 'cp '+pdb_dir+name_array[0]+'.atom'+ ' ' +out_dir
    cmd_2 = 'cp ' + pdb_dir + name_array[1] +'.atom'+ ' ' + out_dir
    print(cmd_1)
    print(cmd_2)
    os.system(cmd_1)
    os.system(cmd_2)
    print(missing_count)