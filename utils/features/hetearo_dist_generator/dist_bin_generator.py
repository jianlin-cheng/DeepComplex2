import copy
import os

import numpy as np


## outputs level L*L where L is the combined length

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


def get_bin_val(_input):
    input_val = float(_input)
    bin_array = np.zeros(37)
    if float(input_val) > 6 or float(input_val) == 0:
        bin_array[0] = 1
    if float(input_val) == 0:
        return bin_array
    if float(input_val) >= 19.5:
        bin_array[36] = 1
        return bin_array
    if float(input_val) < 2.5:
        bin_array[1] = 1
    elif float(input_val) == 6.0:
        #forciing it to be in the 8th bin
        bin_array[8] = 1
    else:
        #
        temp_index = int((input_val / 0.5)) - 4 + 2 - 1
        bin_array[temp_index] = 1
    return bin_array


def label_formatter(_arr, _len):
    formatted_array = np.zeros((_len, _len, 37))
    for val in _arr:
        bin_dist = get_bin_val(val[3])
        formatted_array[val[0] - 1][int(val[1] - 1)] = bin_dist
        formatted_array[int(val[1] - 1)][val[0] - 1] = bin_dist
    return formatted_array


def file_reader(_input):
    content_arry = []
    f = open(_input, "r")
    if f.mode == 'r':
        content_arry = f.read().splitlines()
        f.close()
    return content_arry


fasta_dict = loadFastaDictionary('/home/rajroy/Downloads/fasta_dictionary.txt')
file_arr = []


# input_files = '/home/rajroy/test_dist_dir/all_training_protein_list.txt'
input_files = '/home/rajroy/het_30_tr_roseeta_data_dncon2_divi_feature/train_list_0116221/all_training_protein_list.txt'
list_files = file_reader(input_files)
# output_dir = "/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/rajroy/het_30_bin_200/"
output_dir = "/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/rajroy/NEW_DIST_FILES/dist_bins/"
# dist_dir= "/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/rajroy/het_30_dist/"
dist_dir= "/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/rajroy/NEW_DIST_FILES/Dist_files/"
print(list_files)
missing_counter = 0
for file in list_files:
    name = []
    if '__' in file:
        # temp=files.replace('__','_')
        # name.append(temp.split('_')[1])
        # name.append(temp.split('_')[0])
        name_arr = file.replace('__', '_').split('_')
        name = name_arr[1] + '_' + name_arr[0]
        order_a = 0
        order_b = len(fasta_dict.get(name_arr[1]))

    else:
        # name =files.split('_')
        name_arr = file.split('_')
        name = name_arr[0] + '_' + name_arr[1]
        order_a = 0
        order_b = len(fasta_dict.get(name_arr[0]))

    len_a = len(fasta_dict.get(name_arr[0]))
    len_b = len(fasta_dict.get(name_arr[1]))
    total = len_a + len_b
    # print(str(len_a) + ',' + str(len_b) + ',' + str(len_a + len_b))
    true_file_format = dist_dir+name + '_dist_.txt'
    # get length
    final_name = output_dir + name + '.npz'
    print(file)
    if os.path.exists(true_file_format) and not os.path.exists(final_name):
        missing_counter = missing_counter + 1
        dist_array = read_pair_file_into_array(true_file_format)
        filter_array = filter_dist_array(dist_array, order_a, order_b)
        formatted_value = label_formatter(filter_array, total)
        np.savez_compressed(final_name, formatted_value)
print(missing_counter)
