import os
import sys
#can be used as both len list and list maker

#hetero all len list_maker
def write2file(file, contents):
    with open(file, "w") as f:
        f.writelines(contents)

def single_fasta_file_reader(_seq_file):
    file = open(_seq_file, "r")
    output_array = []
    if file.mode == 'r':
        output_array = file.read().splitlines()
        file.close()
    return output_array




def fasta_file_reader(_seq_file):
    file = open(_seq_file, "r")
    temp = ''
    if file.mode == 'r':
        temp = file.read().splitlines()[0]
        file.close()
    return temp.split(' ')[1], temp.split(' ')[2]


def specific_dir_reader(_input_dir):
    file_names = []
    i = 0
    for root, directories, files in os.walk(_input_dir):
        i = i + 1
        file_names.append(directories)

    return file_names[0]


def specific_filename_reader(_input_dir, _extension):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file:
                file_names.append(file.split(".")[0])
    return file_names


def loadFastaDictionary(dict_file):
    fasta_dict = {}
    #    i=0
    with open(dict_file, "r") as f:
        for line in f:
            #            i+=1
            fasta_dict[line.strip().split(":")[0].strip()] = line.strip().split(":")[1].strip()

    #    print (len(fasta_dict.keys()),i)
    return fasta_dict


fasta_dir = "/home/rajroy/Downloads/fasta_dictionary.txt"
fasta_dict = loadFastaDictionary(fasta_dir)

dir_name = '/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/rajroy/ALIGNMENT_HETERO_30/OUTPUT/LEWIS/a3m_400_more/only_400_a3m_features/'
# dir_name=sys.argv[1]
# name = sys.argv[2]
string_len_list = ''
dir_names = specific_filename_reader(dir_name, ".npz")
counter = 0
lowest_fas= 0
for dir in dir_names:

    name = dir.replace('__','_').split('_')
    name_a = name[0]
    name_b = name[1]
    fasta_a=fasta_dict.get(name_a)
    fasta_b=fasta_dict.get(name_b)
    print(len(fasta_a))
    print(len(fasta_b))
    if '__' in dir:
        final_name =  name_a = name[0]+'__'+ name[1]+'_'+name[2]
    else:
        final_name = name_a = name[0] + '_' + name[1] + '_' + name[2]
    if len(fasta_b)<30 or len(fasta_a)<30:
        counter=counter+1
    else:
        ###len list
        string_len_list=string_len_list+final_name+'\t'+str(len(fasta_a)+len(fasta_b))+"\n"
        ###list
        # string_len_list=string_len_list+final_name+ "\n"
print(counter)
write2file('/home/rajroy/all_list.txt',string_len_list)