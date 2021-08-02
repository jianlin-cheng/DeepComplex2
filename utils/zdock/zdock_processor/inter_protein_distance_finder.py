import os


def specific_dir_reader(_input_dir):
    file_names = []
    i = 0
    for root, directories, files in os.walk(_input_dir):
        i = i + 1
        file_names.append(directories)

    return file_names[0]


_input_pdb = "/home/rajroy/400_zdock_out/val/pdbs/"
list_of_pdb_dir = specific_dir_reader(_input_pdb)
output_dir = '/home/rajroy/400_zdock_out/rr_test/'
pdb_dir = '/home/rajroy/400_zdock_out/seperated_pdbs/'

program= "/home/rajroy/pdb2distance_inter_heavy.py"
for pdb in list_of_pdb_dir:

    names = pdb.split("_")
    file_name_a = pdb_dir + "/" + names[0]+".atom"
    file_name_b = pdb_dir + "/" + names[1]+".atom"
    
    #### sample command
    # ##### python pdb2distance_inter_heavy.py ./1AVYA.atom ./1AVYB.atom 6.0 ./inter_1AVY_AB.txt
    cmd = "python "+program +' '+file_name_a +' '+file_name_b +' 6.0 '+output_dir+pdb+'.rr' +' & '
    print(cmd)
    # os.system(cmd)
