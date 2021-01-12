import copy
import os
import random
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

def specific_filename_reader_with_filters(_input_dir, _extension,_feature_dir):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file:
                if os.path.getsize(_input_dir+'/'+file)>2000:
                    if os.path.isfile(_feature_dir+'/'+file.replace(".a3m",".npz")):
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


def training_list_maker_(_array,_code=0 , _dir=""):
    #all_training_protein_length
    #all_training_protein_list
    #test_list
    #validation_list
    #train_list
    string_len_list=""
    counter=0
    dir_names= copy.deepcopy(_array)
    for dir in dir_names:

        name = dir.replace('__', '_').split('_')
        name_a = name[0]
        name_b = name[1]
        fasta_a = fasta_dict.get(name_a)
        fasta_b = fasta_dict.get(name_b)
        print(len(fasta_a))
        print(len(fasta_b))
        if '__' in dir:
            final_name = name_a = name[0] + '__' + name[1] + '_' + name[2]
        else:
            final_name = name_a = name[0] + '_' + name[1] + '_' + name[2]
        if len(fasta_b) < 30 or len(fasta_a) < 30 or (len(fasta_b) + len(fasta_a)) > 400:
            counter = counter + 1
        else:
            ###len list
            if _code == 0:
                string_len_list=string_len_list+final_name+ "\n"
            else:
                string_len_list = string_len_list + final_name + '\t' + str(len(fasta_a) + len(fasta_b)) + "\n"
            ###list

    if _code == 0:
        array = string_len_list.splitlines()
        test_list_len = len(array)/10
        total_list_len = len(array)
        write2file(_dir + 'all_training_protein_list.txt', string_len_list)

        counter = 0
        str_rows = ""
        train_index= int(total_list_len-(2*test_list_len))
        test_index=int(total_list_len-test_list_len)
        for val in array[0:train_index]:
            str_rows=str_rows+val+"\n"
        write2file(_dir + 'train_list.txt',str_rows)
        str_rows=""
        for val in array[train_index:test_index]:
            str_rows=str_rows+val+"\n"
        write2file(_dir + 'validation_list.txt',str_rows)
        str_rows = ""
        for val in array[test_index:len(array)]:
            str_rows = str_rows + val + "\n"
        write2file(_dir + 'test_list.txt', str_rows)
    else:
        write2file(_dir + 'all_training_protein_length.txt', string_len_list)

    return string_len_list
make_train_set=  1
fasta_dir = "/home/rajroy/Downloads/fasta_dictionary.txt"
fasta_dict = loadFastaDictionary(fasta_dir)
a3m_200_dir='/media/rajroy/My Passport/Alignment_Hetero/200/200_a3m_rr_all/'
a3m_400_dir="/media/rajroy/My Passport/Alignment_Hetero/400/400_a3m_rr_all/"
a3m_500_dir="/media/rajroy/My Passport/Alignment_Hetero/500/500_a3m_rr_all/"

feature_dir_200="/media/rajroy/My Passport/Alignment_Hetero/200/200_features_all/"
feature_dir_400="/media/rajroy/My Passport/Alignment_Hetero/400/400_features_all/"
feature_dir_500="/media/rajroy/My Passport/Alignment_Hetero/500/500_feature_all/"


# dir_name = '/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/rajroy/ALIGNMENT_HETERO_30/OUTPUT/LEWIS/a3m_400_more/only_400_a3m_features/'
# dir_name=sys.argv[1]
# name = sys.argv[2]

# dir_names = specific_filename_reader(dir_name, ".npz")
a3m_200=specific_filename_reader_with_filters(a3m_200_dir,'a3m',feature_dir_200)
a3m_400=specific_filename_reader_with_filters(a3m_400_dir,'a3m',feature_dir_400)
a3m_500=specific_filename_reader_with_filters(a3m_500_dir,'a3m',feature_dir_500)
dir_names=[]
dir_names.extend(a3m_200)
dir_names.extend(a3m_400)
dir_names.extend(a3m_500)
counter = 0
lowest_fas= 0
#removing if any duplicates

dir_names = list(dict.fromkeys(dir_names))
#making sure order is random
for i in range(0,1000):
    random.shuffle(dir_names)

#
# training_list_maker_()
# print(counter)
# write2file('/home/rajroy/all_list.txt',string_len_list)

NAME="400_run"
if make_train_set == 1:

    output_dir=  "/home/rajroy/"
    target_dir_name= output_dir+NAME+'/'
    os.system("mkdir -p "+output_dir+NAME)
    training_list_maker_(_array=dir_names,_code=0,_dir=target_dir_name)
    training_list_maker_(_array=dir_names, _code=1, _dir=target_dir_name)
    # print(counter)
    # write2file(target_dir_name+'all_training_protein_list.txt', target_dir_name,0)