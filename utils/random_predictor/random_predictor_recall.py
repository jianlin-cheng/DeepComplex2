import copy
import os
#this is not used and might be deleted in the furue
###### this needs to be fixed
# ### Recall formula is wrong
import numpy as np

from utils.evalutaion.relaxed_cmaps import make_relax


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
    if true_positive > 0:
        return 100 * true_positive / number_contacts
    else:
        return 0


THRESHOLD = 0.15
fasta_dict = loadFastaDictionary('/home/rajroy/Downloads/experiment_batch/fasta_dictionary.txt')
test_file = '/home/rajroy/het_30_dncon2_model_tr_roseeta_v3_new/training_list_het_121220/validation_list.txt'
# recall just extract how many
cmap_dir = "/home/rajroy/predict_cmap_200_hetero/"
test_file_name = file_reader(test_file)
val_array = []
all_threshold_values = []
for n in range(5, 100, 5):
    THRESHOLD = n / 100
    for file in test_file_name:
        temp_relax_avalue_array = []

        predict_cmap = cmap_dir + file + '_.rr.npy.txt'
        real_cmap = cmap_dir + 'Y-' + file + '.txt.npy.txt'

        name_arrr = os.path.basename(predict_cmap).split('_')

        len_a = len(fasta_dict.get(name_arrr[0]))
        len_b = len(fasta_dict.get(name_arrr[1]))
        total = len_b + len_a

        real_arr = getY(real_cmap)
        random_cmap = np.random.random((total, total))
        new_pred_map = fix_pred_map(random_cmap, len_a, len_b)
        new_pred_map[new_pred_map > THRESHOLD] = 1
        temp_relax_avalue_array.append(calcuate_recall(new_pred_map, real_arr))
        real_arr_1 = make_relax(real_arr, 1)
        temp_relax_avalue_array.append(calcuate_recall(new_pred_map, real_arr_1))
        real_arr_2 = make_relax(real_arr, 2)
        temp_relax_avalue_array.append(calcuate_recall(new_pred_map, real_arr_2))
        val_array.append(temp_relax_avalue_array)

    sum_relax_0 = 0
    sum_relax_1 = 0
    sum_relax_2 = 0
    for values in val_array:
        sum_relax_0 += values[0]
        sum_relax_1 += values[1]
        sum_relax_2 += values[2]


    temp_all_threshold_values = [THRESHOLD, sum_relax_0 / len(val_array), sum_relax_1 / len(val_array),
                                 sum_relax_2 / len(val_array)]
    all_threshold_values.append(temp_all_threshold_values)


print('THRESHOLD' + '\t\t' + '  RELAX 0  ' + '\t\t' + '  RELAX 1  ' + '\t\t' + '  RELAX 2  ' + '\t\t')
for values in all_threshold_values:
    print(str(values[0]) + '\t\t' + str(values[1]) + '\t\t' + str(values[2]) + '\t\t' + str(values[3]) + '\t\t')
