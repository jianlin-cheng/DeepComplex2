import numpy as np


def loadFastaDictionary(dict_file):
    fasta_dict = {}
    with open(dict_file, "r") as f:
        for line in f:
            fasta_dict[line.strip().split(":")[0].strip()] = line.strip().split(":")[1].strip()
    return fasta_dict


# get value
# append to both

def binding_only_file_reader(_seq_file):
    my_dictionary = {}
    with open(_seq_file) as infile:
        number_counter = 0
        for line in infile:
            # print(line)
            # getting the binding file

            temp_array = line.split("\t")
            list_a = []
            list_b = []
            # fasta name
            fasta_name_a = temp_array[0]
            fasta_name_b = temp_array[1]

            if my_dictionary.get(fasta_name_a) != None:
                list_a = my_dictionary[fasta_name_a]

            if my_dictionary.get(fasta_name_b) != None:
                list_b = my_dictionary[fasta_name_b]

                # get list
            list_a.append(fasta_name_b)
            list_b.append(fasta_name_a)
            # append list

            #
            my_dictionary[fasta_name_a] = list(set(list_a))
            my_dictionary[fasta_name_b] = list(set(list_b))
            number_counter = number_counter + 1
            if number_counter % 1000000 == 0:
                print('number_counter ' + str(number_counter))
    np.save('STRINGS_BINDING_DB_DICT.npy', my_dictionary)


binding_only_file_reader('/exports/store1/raj/future_task_new_db/binding_file.txt')
# np.load('/home/rajroy/STRINGS_DB_DICT.npy', allow_pickle='TRUE').item()
