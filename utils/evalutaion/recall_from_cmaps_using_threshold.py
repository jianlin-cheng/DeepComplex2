import copy
import os

import numpy as np

from utils.evalutaion.relaxed_cmaps import make_relax

# # NOTE
## This basically takes 2 contact maps as input and uses them to find the recall
## converts the real one into relax map of 0,1,2
#######  Threshold was used
###### this needs to be fixed
# ### Recall formula is wrong

def file_reader(_input):
    content_arry = []
    f = open(_input, "r")
    if f.mode == 'r':
        content_arry = f.read().splitlines()
        f.close()
    return content_arry


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
            output_rr[j_counter][i_counter] = (output_rr[i_counter][j_counter] + output_rr[j_counter][i_counter]) / 2.0
            output_rr[i_counter][j_counter] = 0
    return output_rr


def getY(true_file):
    # calcualte the length of the protein (the first feature)
    input_array = []
    file = open(true_file, "r")
    if file.mode == 'r':
        input_array = file.read().splitlines()
        file.close()
    inter_array = []

    for values in input_array:
        inter_array.append(values.strip().split(' '))
    # since its a sequare matrix coz I made it that way
    L = len(inter_array)
    relax_0_array = np.asfarray(inter_array, float)
    # print(inter_array)

    return relax_0_array


def loadFastaDictionary(dict_file):
    fasta_dict = {}
    with open(dict_file, "r") as f:
        for line in f:
            fasta_dict[line.strip().split(":")[0].strip()] = line.strip().split(":")[1].strip()
    return fasta_dict


def calcuate_recall(_pred, _true):
    pred_eval = np.where(_pred == 1.0)
    number_contacts = len(pred_eval[0])
    true_positive = 0
    for val in range(0, number_contacts):
        index_X = pred_eval[0][val]
        index_Y = pred_eval[1][val]
        if _true[index_X][index_Y] == 1.0:
            true_positive = true_positive + 1
    if true_positive >0:
        return 100 * true_positive / number_contacts
    else:
        return 0

THRESHOLD = 0.15
fasta_dict = loadFastaDictionary('/home/rajroy/Downloads/experiment_batch/fasta_dictionary.txt')
test_file = '/home/rajroy/het_30_dncon2_model_tr_roseeta_v3_new/training_list_het_121220/validation_list.txt'
# recall just extract how many
cmap_dir = "/home/rajroy/predict_cmap_200_hetero/"
test_file_name = file_reader(test_file)
val_array= []
all_threshold_values = []
for n in range(5,100,5):
    THRESHOLD = n/100
    for file in test_file_name:
        temp_relax_avalue_array=[]
        predict_cmap = cmap_dir + file + '_.rr.npy.txt'
        # /home/rajroy/predict_cmap_200_hetero/1H3OB_1H3OC_3305_.rr.npy.txt'
        real_cmap = cmap_dir + 'Y-' + file + '.txt.npy.txt'
        # /home/rajroy/predict_cmap_200_hetero/Y-1H3OB_1H3OC_3305.txt.npy.txt'
        pred_arr = getY(predict_cmap)
        # processed cmap

        name_arrr = os.path.basename(predict_cmap).split('_')
        len_a = len(fasta_dict.get(name_arrr[0]))
        len_b = len(fasta_dict.get(name_arrr[1]))
        real_arr = getY(real_cmap)
        new_pred_map = fix_pred_map(pred_arr, len_a, len_b)
        new_pred_map[new_pred_map > THRESHOLD] = 1
        temp_relax_avalue_array.append(calcuate_recall(new_pred_map, real_arr))
        real_arr_1 = make_relax(real_arr, 1)
        temp_relax_avalue_array.append(calcuate_recall(new_pred_map, real_arr_1))
        real_arr_2 = make_relax(real_arr, 2)
        temp_relax_avalue_array.append(calcuate_recall(new_pred_map, real_arr_2))
        val_array.append(temp_relax_avalue_array)

    sum_relax_0=0
    sum_relax_1=0
    sum_relax_2=0
    for values in val_array:
        sum_relax_0+=values[0]
        sum_relax_1 += values[1]
        sum_relax_2 += values[2]
        # print('relax 0 '+str(values[0]) +' relax 1 '+str(values[1])+' relax 2 '+str(values[2])+'\n')

    temp_all_threshold_values=  [THRESHOLD,sum_relax_0/len(val_array),sum_relax_1/len(val_array),sum_relax_2/len(val_array)]
    all_threshold_values.append(temp_all_threshold_values)
    # print('AVERAGE'+'\n')
    # print(str(sum_relax_0/len(val_array))+' '+str(sum_relax_1/len(val_array))+' '+str(sum_relax_2/len(val_array)))

print('THRESHOLD'+'\t\t'+'  RELAX 0  '+'\t\t'+'  RELAX 1  '+'\t\t'+'  RELAX 2  '+'\t\t')
for values in all_threshold_values:
    print(str(values[0])+'\t\t'+str(values[1])+'\t\t'+str(values[2])+'\t\t'+str(values[3])+'\t\t')
