import copy
import os

import numpy as np
## outputs level L*L where L is the combined length


def read_cmap_file_into_array(_input):
    output_array = []
    file = open(_input, "r")
    if file.mode == 'r':
        temp_array = file.read().splitlines()
        file.close()

    for values in temp_array:
        dist_row = []
        temp = values.split(" ")
        output_array.append(temp)
    return output_array



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
def specific_filename_reader(_input_dir, _extension):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file:
                file_names.append(file.split(".")[0])
    return file_names


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

def fix_pred_map(_inputs, _len_a, _len_b):
    total = _len_a + _len_b
    output_rr = copy.deepcopy(_inputs)
    for j_counter in range(0, _len_a):
        for i_counter in range(0, _len_a):
            output_rr[j_counter][i_counter] = 0
            output_rr[i_counter][j_counter] = 0

    for j_counter in range(_len_a, total):
        for i_counter in range(_len_a, total):
            output_rr[j_counter][i_counter] = 0
            output_rr[i_counter][j_counter] = 0

    for j_counter in range(_len_a, total):
        for i_counter in range(0, _len_a):
            output_rr[j_counter][i_counter] = (float(output_rr[i_counter][j_counter]) + float(output_rr[j_counter][i_counter])) / 2.0
            output_rr[i_counter][j_counter] = 0

    return output_rr
def get_cmap_string(_array):
    out_str =""
    for j_counter in range(0, _array.shape[0]):
        temp_str=""
        for k_counter in range(0, _array.shape[1]):
            temp_str = temp_str+ str(_array[j_counter][k_counter])+ " "
        out_str = out_str+temp_str+"\n"
    return out_str

def get_rr_string(_array):
    out_str =""
    for j_counter in range(0, _array.shape[0]):

        for k_counter in range(0, _array.shape[1]):
            line = str(j_counter+1)+" "+str(k_counter+1)+" 0 6 "+str(_array[j_counter][k_counter])
            out_str = out_str+ line+"\n"

    return out_str


def get_filtered_cmap(_array , _len_a,_len_b):
    # out_arry = np.zeros((_len_b,_len_a))
    # for j_counter in range (0,_len_b):
    #     for k_counter in range (0,_len_a):
    #         out_arry [j_counter][k_counter]=_array[j_counter+_len_a][k_counter]
    out_arry = np.zeros((_len_a,_len_b))
    for j_counter in range (0,_len_a):
        for k_counter in range (0,_len_b):
            # if k_counter == 1272:
            #     print("bug")
            out_arry [j_counter][k_counter]=_array[j_counter][k_counter+len_a]
            # print(str(j_counter)+" "+str(k_counter))
    return  out_arry

fasta_dict = loadFastaDictionary("/home/rajroy/msa_hetero_algae/fasta_dict_algae.txt")
out_dir = "/home/rajroy/msa_hetero_algae/cmap/"
rr_dir = "/home/rajroy/algae_1300/"
fasta_lists =specific_filename_reader(_input_dir="/home/rajroy/msa_hetero_algae/concat_fastas/",_extension="fasta")
for file in fasta_lists:
    temp_name=     rr_dir +file+"_.rr.npy.txt"
    final_cmap_name = out_dir +file+"_.cmap_het"
    final_rr_name = out_dir +file+"_.rr_het"
    print(temp_name)
    name_arr = file.split('_')
    # order_a = 0
    # order_b = len(fasta_dict.get(name_arr[0]))
    len_a = len(fasta_dict.get(name_arr[0]))
    len_b = len(fasta_dict.get(name_arr[1]))
    # if not os.path.isfile(final_cmap_name):
    pred_arr = read_cmap_file_into_array(temp_name)
    # new_pred_map = fix_pred_map(pred_arr, len_a, len_b)
    new_pred_map = np.transpose(fix_pred_map(pred_arr, len_a, len_b))
    final_cmap  = get_filtered_cmap(new_pred_map,len_a, len_b)
    print("here")
    contents = get_cmap_string(final_cmap)
        # write2file(final_cmap_name, contents)
        #
        # contents_rr = get_rr_string(final_cmap)
        # write2file(final_rr_name, contents_rr)

