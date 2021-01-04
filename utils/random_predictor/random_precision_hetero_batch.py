import copy
import math
import os

import numpy as np

from utils.evalutaion.relaxed_cmaps import make_relax
### calculates precision with a random contact map

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


def calculateEvaluationStats(_pred_cmap, _true_cmap, L):
    pred_cmap=copy.deepcopy(_pred_cmap)
    true_cmap=copy.deepcopy(_true_cmap)
    prec_T5, prec_T10, prec_T20, prec_T30, prec_T50, prec_L30, prec_L20, prec_L10, prec_L5, prec_L, prec_2L, con_num = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    max_Top = int((2*L ) + 0.5)
    pred_eval = np.where(_true_cmap == 1.0)
    number_contacts = len(pred_eval[0])
    if 50 > max_Top: max_Top = 50

    for i in range(1, max_Top + 1):
        (x, y) = np.unravel_index(np.argmax(pred_cmap, axis=None), pred_cmap.shape)
        pred_cmap[x][y] = 0
        if true_cmap[x][y] == 1:
            con_num += 1
        if i == 5:
            prec_T5 = con_num * 20
            if prec_T5 > 100: prec_T5 = 100
            print("L=", L, "Val=", 5, "Con_num=", con_num)
        if i == 10:
            prec_T10 = con_num * 10
            if prec_T10 > 100: prec_T10 = 100
            print("L=", L, "Val=", 10, "Con_num=", con_num)
        if i == 20:
            prec_T20 = con_num * 5
            if prec_T20 > 100: prec_T20 = 100
            print("L=", L, "Val=", 20, "Con_num=", con_num)
        if i == 30:
            prec_T30 = con_num * 100 / 30
            if prec_T30 > 100: prec_T30 = 100
            print("L=", L, "Val=", 30, "Con_num=", con_num)
        if i == 50:
            prec_T50 = con_num * 2
            if prec_T50 > 100: prec_T50 = 100
            print("L=", L, "Val=", 50, "Con_num=", con_num)
        if i == int((L / 30) + 0.5):
            prec_L30 = con_num * 100 / i
            if prec_L30 > 100: prec_L30 = 100
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((L / 20) + 0.5):
            prec_L20 = con_num * 100 / i
            if prec_L20 > 100: prec_L20 = 100
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((L / 10) + 0.5):
            prec_L10 = con_num * 100 / i
            if prec_L10 > 100: prec_L10 = 100
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((L / 5) + 0.5):
            prec_L5 = con_num * 100 / i
            if prec_L5 > 100: prec_L5 = 100
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((L) + 0.5):
            prec_L = con_num * 100 / i
            if prec_L > 100: prec_L = 100
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((2 * L) + 0.5):
            prec_2L = con_num * 100 / i
            if prec_2L > 100: prec_2L = 100
            print("L=", L, "Val=", i, "Con_num=", con_num)

    return [prec_T5/number_contacts, prec_T10/number_contacts, prec_T20/number_contacts, prec_T30/number_contacts, prec_T50/number_contacts, prec_L30/number_contacts, prec_L20/number_contacts, prec_L10/number_contacts, prec_L5/number_contacts, prec_L/number_contacts, prec_2L/number_contacts]



def get_evaluation_result(_arr,_relax,_SAMPLE_SIZE):
    sum_prec_T5, sum_prec_T10, sum_prec_T20, sum_prec_T30, sum_prec_T50, sum_prec_L30, sum_prec_L20, sum_prec_L10, sum_prec_L5, sum_prec_L, sum_prec_2L = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    SAMPLE_SIZE = len(test_file_name)
    for values in _arr:
        sum_prec_T5 = sum_prec_T5 + values[0]
        sum_prec_T10 = sum_prec_T10 + values[1]
        sum_prec_T20 = sum_prec_T20 + values[2]
        sum_prec_T30 = sum_prec_T30 + values[3]
        sum_prec_T50 = sum_prec_T50 + values[4]
        sum_prec_L30 = sum_prec_L30 + values[5]
        sum_prec_L20 = sum_prec_L20 + values[6]
        sum_prec_L10 = sum_prec_L10 + values[7]
        sum_prec_L5 = sum_prec_L5 + values[8]
        sum_prec_L = sum_prec_L + values[9]
        sum_prec_2L = sum_prec_2L + values[10]
    print(    str(_relax) + '\t\t\t' + str(sum_prec_T5 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_T10 / _SAMPLE_SIZE)[
                                                                               0:5] + '\t\t\t' + str(
            sum_prec_T20 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_T30 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
            sum_prec_T50 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_L30 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
            sum_prec_L20 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_L10 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
            sum_prec_L5 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_L / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
            sum_prec_2L / _SAMPLE_SIZE)[0:5] )


relax_0 = []
relax_1 = []
relax_2 = []
fasta_dict = loadFastaDictionary('/home/rajroy/Downloads/experiment_batch/fasta_dictionary.txt')
test_file = '/home/rajroy/het_30_dncon2_model_tr_roseeta_v3_new/training_list_het_121220/test_list.txt'
# recall just extract how many
cmap_dir = "/home/rajroy/predict_cmap_200_hetero_test/"
test_file_name = file_reader(test_file)
val_array = []
all_threshold_values = []
# for n in range(5,100,5):
# THRESHOLD = n/100
SAMPLE_SIZE=0
for file in test_file_name:
    predict_cmap = cmap_dir + file + '_.rr.npy.txt'
    if os.path.isfile(predict_cmap):
        SAMPLE_SIZE=SAMPLE_SIZE+1
        # if '__' in file:
        #     print(file)
        temp= file
        if '__' in file:
            temp = file.replace('__', '_')

        name_arrr = temp.split('_')


        if '__' in os.path.basename(predict_cmap):
            len_b = len(fasta_dict.get(name_arrr[0]))
            len_a = len(fasta_dict.get(name_arrr[1]))
        else:
            len_a = len(fasta_dict.get(name_arrr[0]))
            len_b = len(fasta_dict.get(name_arrr[1]))
    # /home/rajroy/predict_cmap_200_hetero/1H3OB_1H3OC_3305_.rr.npy.txt'
        real_cmap = cmap_dir + 'Y-' + file + '.txt.npy.txt'
        # /home/rajroy/predict_cmap_200_hetero/Y-1H3OB_1H3OC_3305.txt.npy.txt'
        # pred_arr = getY(predict_cmap)
        # processed cmap

        # name_arrr = os.path.basename(predict_cmap).split('_')
        # len_a = len(fasta_dict.get(name_arrr[0]))
        # len_b = len(fasta_dict.get(name_arrr[1]))
        real_arr = getY(real_cmap)
        # total =math.sqrt(len_a*len_b)
        total = len_a + len_b
        random_cmap = np.random.random((total, total))
        new_pred_map = fix_pred_map(random_cmap, len_a, len_b)
        relax_0.append(calculateEvaluationStats(new_pred_map, real_arr, total))

        real_arr_1 = make_relax(real_arr, 1)
        relax_1.append(calculateEvaluationStats(new_pred_map, real_arr_1, total))

        real_arr_2 = make_relax(real_arr, 2)
        relax_2.append(calculateEvaluationStats(new_pred_map, real_arr_2, total))

print(
    'RELAX ' + '\t\t\t' + 'TOP-5' + '\t\t\t' + 'TOP-10' + '\t\t\t' + 'TOP-20' + '\t\t\t' + 'TOP-30' + '\t\t\t' + 'TOP-50' + '\t\t\t' + 'L/30' + '\t\t\t' + 'L/20' + '\t\t\t' + 'L/10' +
    '\t\t\t' + 'L/5' + '\t\t\t' + 'L/' + '\t\t\t' + '2L/')

get_evaluation_result(relax_0,0,len(test_file_name))
get_evaluation_result(relax_1,1,len(test_file_name))
get_evaluation_result(relax_2,2,len(test_file_name))
