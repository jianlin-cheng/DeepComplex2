import copy
import math
import os

import numpy as np

from utils.evalutaion.relaxed_cmaps import make_relax
# # NOTE
## This basically takes 2 contact maps as input and uses them to find the recall
## converts the real one into relax map of 0,1,2
### *For Relax 0,1 and 2 only the true contact map was taken into consideration not the relaxed maps to identify the false negative

####### Recall = true_positive /(true_positve +false negative)
#####  E.g For calculating recall of top-5 for a protein of 10 true contacts, if we predict 3 contacts successfully, then
#### True positive  = 3
#### False negative = 10-3 = 7 , as 7 true contacts were not predicted correctly
#### Recall = 3/ (3 +7 ) =0.3
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


def calculateEvaluationStats(_pred_cmap, _true_cmap, L,_number_contacts):
    pred_cmap=copy.deepcopy(_pred_cmap)
    true_cmap=copy.deepcopy(_true_cmap)
    prec_T5, prec_T10, prec_T20, prec_T30, prec_T50, prec_L30, prec_L20, prec_L10, prec_L5, prec_L, prec_2L, con_num = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    max_Top = int((2*L ) + 0.5)
    if 50 > max_Top: max_Top = 50

    for i in range(1, max_Top + 1):
        (x, y) = np.unravel_index(np.argmax(pred_cmap, axis=None), pred_cmap.shape)
        pred_cmap[x][y] = 0
        if true_cmap[x][y] == 1:
            con_num += 1

        if i == 5:
            prec_T5 = con_num/_number_contacts
            if prec_T5  > 1.0:
                prec_T5 = 1.0
            print("L=", L, "Val=", 5, "Con_num=", con_num)
        if i == 10:
            prec_T10 = con_num/_number_contacts
            if prec_T10  > 1.0: prec_T10 = 1.0
            print("L=", L, "Val=", 10, "Con_num=", con_num)
        if i == 20:
            prec_T20 = con_num/_number_contacts
            if prec_T20  > 1.0: prec_T20 = 1.0
            print("L=", L, "Val=", 20, "Con_num=", con_num)
        if i == 30:
            prec_T30 = con_num/_number_contacts
            if prec_T30  > 1.0: prec_T30 = 1.0
            print("L=", L, "Val=", 30, "Con_num=", con_num)
        if i == 50:
            prec_T50 = con_num/_number_contacts
            if prec_T50  > 1.0: prec_T50 = 1.0
            print("L=", L, "Val=", 50, "Con_num=", con_num)
        if i == int((L / 30) + 0.5):
            prec_L30 = con_num/_number_contacts
            if prec_L30  > 1.0: prec_L30 = 1.0
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((L / 20) + 0.5):
            prec_L20 = con_num/_number_contacts
            if prec_L20  > 1.0: prec_L20 = 1.0
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((L / 10) + 0.5):
            prec_L10 = con_num/_number_contacts
            if prec_L10  > 1.0: prec_L10 = 1.0
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((L / 5) + 0.5):
            prec_L5 = con_num/_number_contacts
            if prec_L5  > 1.0: prec_L5 = 1.0
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((L) + 0.5):
            prec_L = con_num/_number_contacts
            if prec_L  > 1.0: prec_L = 1.0
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((2 * L) + 0.5):
            prec_2L = con_num/_number_contacts
            if prec_2L  > 1.0: prec_2L = 1.0
            print("L=", L, "Val=", i, "Con_num=", con_num)

    # return [prec_T5/_number_contacts, prec_T10/_number_contacts, prec_T20/_number_contacts, prec_T30/_number_contacts, prec_T50/_number_contacts, prec_L30/_number_contacts, prec_L20/_number_contacts, prec_L10/_number_contacts, prec_L5/_number_contacts, prec_L/_number_contacts, prec_2L/_number_contacts]
    return [prec_T5, prec_T10, prec_T20, prec_T30, prec_T50, prec_L30, prec_L20, prec_L10, prec_L5, prec_L, prec_2L]



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
# test_file = '/home/rajroy/het_30_dncon2_model_tr_roseeta_v3_new/training_list_het_121220/test_list.txt'
# test_file = '/home/rajroy/400_run/test_list.txt'
# recall just extract how many
# cmap_dir = "/home/rajroy/hq_200/"
# cmap_dir = "/home/rajroy/400_new_het_24/"
# cmap_dir = "/home/rajroy/cmap_predict_200_farhan/"

test_file = '/home/rajroy/het_30_dncon2_model_tr_roseeta_v3_new_2/training_list_het_400/test_list.txt'
# test_file = '/home/rajroy/het_30_dncon2_model_tr_roseeta_v3_new/training_list_het_400/validation_list.txt'
# cmap_dir = "/home/rajroy/cmap_predict_400_hetero/"
# cmap_dir = "/home/rajroy/400_new_het_24/"
# cmap_dir = "/home/rajroy/predict_cmap_200_hetero_test/"
cmap_dir = "/home/rajroy/400_zdock_out/cmaps/"
test_file_name = file_reader(test_file)
val_array = []
all_threshold_values = []
# for n in range(5,100,5):
# THRESHOLD = n/100
SAMPLE_SIZE=0
for file in test_file_name:
    name_1 = file.split('_')
    predict_cmap = cmap_dir + name_1[0]+'_'+name_1[1] + '.txt'


    # predict_cmap = cmap_dir + name_arrr_1[0]+'_'+name_arrr_1[1] + '.txt'
    if os.path.isfile(predict_cmap):
        SAMPLE_SIZE=SAMPLE_SIZE+1

        # /home/rajroy/predict_cmap_200_hetero/1H3OB_1H3OC_3305_.rr.npy.txt'

        # /home/rajroy/predict_cmap_200_hetero/Y-1H3OB_1H3OC_3305.txt.npy.txt'

        # processed cmap
        name = os.path.basename(predict_cmap).replace('.txt','')
        if '__' in name:
            name=name.replace('__','_')

        name_arrr = name.split('_')
        if '__' in  os.path.basename(predict_cmap):
            len_b = len(fasta_dict.get(name_arrr[0]))
            len_a = len(fasta_dict.get(name_arrr[1]))
        else:
            len_a = len(fasta_dict.get(name_arrr[0]))
            len_b = len(fasta_dict.get(name_arrr[1]))
        total = len_a + len_b

        #removing symmetry
        # real_cmap = cmap_dir + 'Y-' + file + '.txt.npy.txt'
        real_cmap = cmap_dir + 'Y-' + file + '.txt.npy.txt'
        real_arr = fix_pred_map( getY(real_cmap), len_a, len_b)

        #getting number of contacts
        pred_eval = np.where(real_arr == 1.0)
        number_contacts = len(pred_eval[0])


        random_cmap = np.random.random((total, total))
        pred_arr = getY(predict_cmap)
        new_pred_map = fix_pred_map(pred_arr, len_a, len_b)

        relax_0.append(calculateEvaluationStats(new_pred_map, real_arr, total,number_contacts))

        real_arr_1 = make_relax(real_arr, 1)
        relax_1.append(calculateEvaluationStats(new_pred_map, real_arr_1, total,number_contacts))

        real_arr_2 = make_relax(real_arr, 2)
        relax_2.append(calculateEvaluationStats(new_pred_map, real_arr_2, total,number_contacts))

print(
    'RELAX ' + '\t\t\t' + 'TOP-5' + '\t\t\t' + 'TOP-10' + '\t\t\t' + 'TOP-20' + '\t\t\t' + 'TOP-30' + '\t\t\t' + 'TOP-50' + '\t\t\t' + 'L/30' + '\t\t\t' + 'L/20' + '\t\t\t' + 'L/10' +
    '\t\t\t' + 'L/5' + '\t\t\t' + 'L/' + '\t\t\t' + '2L/')

get_evaluation_result(relax_0,0,len(test_file_name))
get_evaluation_result(relax_1,1,len(test_file_name))
get_evaluation_result(relax_2,2,len(test_file_name))
print(SAMPLE_SIZE)