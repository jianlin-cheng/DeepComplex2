import os


def single_fasta_file_reader(_seq_file):
    file = open(_seq_file, "r")
    if file.mode == 'r':
        output_array = file.read().splitlines()
        file.close()

    final_array = []
    for value in output_array:
        temp = []
        val_a = value.replace('.atom', '')
        temp.append(val_a.split(',')[0])
        temp.append(val_a.split(',')[1])
        final_array.append(temp)

    return final_array


def specific_filename_reader(_input_dir, _extension):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file:
                file = file.split(".")[0]
                if not file in file_names:
                    file_names.append(file.split(".")[0])

    return file_names


z_dir = '/home/rajroy/Downloads/zdock3.0.2_linux_x64/'

test_list_file = '/home/rajroy/paired_test_371_500.txt'
test_pdb_file = '/home/rajroy/atom_files_370_500_test/'
test_out_file = '/home/rajroy/out_dock/test/'


test_pdb_list = single_fasta_file_reader(test_list_file)
# test_pdb_array = sorted(test_pdb_list)
# val_list=specific_filename_reader(test_pdb_file,'atom')

counter_test = 0
for t in test_pdb_list:
    name = t[0][0:4]
    both_name = t[0] + '_' + t[1]
    z_file = test_out_file + '/files/' + name + '.out'
    file_a = test_pdb_file + '/' + t[0] + '.atom'
    file_b = test_pdb_file + '/' + t[1] + '.atom'

    z_cmd = z_dir + '/zdock -R ' + file_a + ' -L ' + file_b + ' -N 5 -o ' + z_file
    if not os.path.exists(z_file):
        val_out = os.system(z_cmd)

        os.chdir(z_dir)
        create_cmd = 'perl ' + 'create.pl ' + z_file
        print(create_cmd)

        os.system(create_cmd)
        target_pdb_file = test_out_file + '/' + 'pdbs/' + both_name
        os.system('mkdir -p ' + target_pdb_file)
        os.system('mv *.pdb ' + target_pdb_file)
        counter_test = counter_test + 1
        print(counter_test)
    # run program
    # create it
    # copy to the designated
#
# val_list_file = '/home/rajroy/homo_std_docking/homo_std_val_pairs.txt'
# val_pdb_file = '/home/rajroy/homo_std_docking/atom_files_val/'
# val_out_file = '/home/rajroy/out_dock/val/'
# val_pdb_list = single_fasta_file_reader(val_list_file)
# # test_pdb_array = sorted(val_pdb_list)
# # val_list=specific_filename_reader(test_pdb_file,'atom')
#
# counter_val = 0
# for t in val_pdb_list:
#     name = t[0][0:4]
#     both_name = t[0] + '_' + t[1]
#     z_file = val_out_file + '/files/' + name + '.out'
#     file_a = val_pdb_file + '/' + t[0] + '.atom'
#     file_b = val_pdb_file + '/' + t[1] + '.atom'
#
#     z_cmd = z_dir + '/zdock -R ' + file_a + ' -L ' + file_b + ' -N 5 -o ' + z_file
#     if not os.path.exists(z_file):
#         val_out = os.system(z_cmd)
#
#         os.chdir(z_dir)
#         create_cmd = 'perl ' + 'create.pl ' + z_file
#         print(create_cmd)
#
#         os.system(create_cmd)
#         target_pdb_file = val_out_file + '/' + 'pdbs/' + both_name
#         os.system('mkdir -p ' + target_pdb_file)
#         os.system('mv *.pdb ' + target_pdb_file)
#         counter_val=counter_val+1
#         print(counter_val)