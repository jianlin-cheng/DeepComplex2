import copy
import os
import random

import numpy as np
#simple program to make a random contact map

def filter_dist_array(_dist, i_adder, j_adder):
    original = copy.deepcopy(_dist)
    out_array = []
    for val in original:
        # temp = []
        temp = [int(val[0]) + i_adder, int(val[1]) + j_adder, val[2], val[3], val[4]]
        out_array.append(temp)
    return out_array


def read_pair_file_into_array(_input):
    output_array = []
    file = open(_input, "r")
    if file.mode == 'r':
        temp_array = file.read().splitlines()
        file.close()

    for values in temp_array:
        dist_row = []
        temp = values.split(" ")
        if len(temp) == 5:
            dist_row.append(temp[0])
            dist_row.append(temp[1])
            dist_row.append(temp[2])
            dist_row.append(temp[3])
            dist_row.append(temp[4])
            output_array.append(dist_row)
    return output_array


def specific_filename_reader(_input_dir, _extension):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file and "start" not in file:
                file_names.append(file)
    return file_names


def len_reader(_input_file):
    true_distance_file = open(_input_file, "r")

    if true_distance_file.mode == 'r':
        output_array = true_distance_file.read()
        true_distance_file.close()

    return len(output_array.splitlines()[0])


def loadFastaDictionary(dict_file):
    fasta_dict = {}
    #    i=0
    with open(dict_file, "r") as f:
        for line in f:
            fasta_dict[line.strip().split(":")[0].strip()] = line.strip().split(":")[1].strip()
    return fasta_dict


def write2file(file, contents):
    with open(file, "w") as f:
        f.writelines(contents)
    f.close()


def file_array_return(_input_dir):
    output_array = []
    true_distance_file = open(_input_dir, "r")

    if true_distance_file.mode == 'r':
        output_array = true_distance_file.read().splitlines()
        true_distance_file.close()

    return output_array


def process_file(_rr_array, _len_a):
    filered_array = []
    # add to i or j
    for val in _rr_array:
        i = int(val[0])
        j = int(val[1])

    return None


def specific_dir_reader(_input_dir):
    file_names = []
    i = 0
    for root, directories, files in os.walk(_input_dir):
        i = i + 1
        file_names.append(directories)

    return file_names[0]


def label_formatter(_arr, _len):
    formatted_array = np.zeros((_len, _len))
    for val in _arr:
        formatted_array[val[0] - 1][int(val[1] - 1)] = 1
        formatted_array[int(val[1] - 1)][val[0] - 1] = 1
    out_string = ''
    for values in formatted_array:
        one_row = ''
        for d in values:
            one_row = one_row + str(d) + str(' ')
        out_string = out_string + one_row + '\n'
    return out_string


fasta_dict = loadFastaDictionary('/home/rajroy/Downloads/fasta_dictionary.txt')

fasta_dir = '/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/rajroy/back_up/het30_tr_roseeta_training_data/het30_splitted_fasta/'
file = fasta_list = specific_dir_reader(fasta_dir)[0]

file_arr = []

total = 0
count = 0
# geta  fastas
# read it from fasta Dicta
# find the order
# make empty numpy
# rnadom values

# for file in fasta_list:
print(file)

# name_arr = file.split('_')
name = ''
order_a = 0
order_b = 0

if '__' in file:

    name_arr = file.replace('__', '_').split('_')
    name = name_arr[1] + '_' + name_arr[0]

    order_b = len(fasta_dict.get(name_arr[1]))

else:
    name_arr = file.split('_')
    name = name_arr[0] + '_' + name_arr[1]

len_a = len(fasta_dict.get(name_arr[0]))
len_b = len(fasta_dict.get(name_arr[1]))
total = len_a + len_b
print(str(len_a) + ',' + str(len_b) + ',' + str(len_a + len_b))
random_cmap = np.random.random((total, total))
print(random_cmap)
