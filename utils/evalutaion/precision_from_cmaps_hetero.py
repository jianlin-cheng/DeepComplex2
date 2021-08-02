import copy
import csv
import math
import os

import numpy as np

from utils.evalutaion.relaxed_cmaps import make_relax

# # NOTE
## This basically takes 2 contact maps as input and uses them to find the precision
## converts the real one into relax map of 0,1,2
##compares them and usese out legacy code to calculate the values

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


def calculateEvaluationStats(_pred_cmap, _true_cmap, L,_name):
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

    return [prec_T5, prec_T10, prec_T20, prec_T30, prec_T50, prec_L30, prec_L20, prec_L10, prec_L5, prec_L, prec_2L,_name]



def get_evaluation_result(_arr,_relax,_SAMPLE_SIZE):
    sum_prec_T5, sum_prec_T10, sum_prec_T20, sum_prec_T30, sum_prec_T50, sum_prec_L30, sum_prec_L20, sum_prec_L10, sum_prec_L5, sum_prec_L, sum_prec_2L = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    SAMPLE_SIZE = len(test_file_name)
    some_array= []
    fasta_dict = loadFastaDictionary('/home/rajroy/Downloads/experiment_batch/fasta_dictionary.txt')
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

        name_arrr=values[11].replace("__","_").split('_')
        if fasta_dict.get(name_arrr[0]) == None:
            print('her')
        len_a = len(fasta_dict.get(name_arrr[0]))
        len_b = len(fasta_dict.get(name_arrr[1]))
        total=len_a+len_b
        ratio = len_a/len_b
        if ratio >1 :
            ratio= 1 /ratio



        some_array.append([values[11],str(values[0])[0:5],str(values[1])[0:5],str(values[2])[0:5],str(values[3])[0:5],str(values[4])[0:5],
                           str(values[5])[0:5],str(values[6])[0:5],str(values[7])[0:5],str(values[8])[0:5],str(values[9])[0:5],str(values[10])[0:5],str(total),str(ratio)[0:3]])

    print(    str(_relax) + '\t\t\t' + str(sum_prec_T5 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_T10 / _SAMPLE_SIZE)[
                                                                               0:5] + '\t\t\t' + str(
            sum_prec_T20 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_T30 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
            sum_prec_T50 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_L30 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
            sum_prec_L20 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_L10 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
            sum_prec_L5 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_L / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
            sum_prec_2L / _SAMPLE_SIZE)[0:5] )

    return some_array
relax_0 = []
relax_1 = []
relax_2 = []
fasta_dict = loadFastaDictionary('/home/rajroy/Downloads/experiment_batch/fasta_dictionary.txt')
# test_file = '//home/rajroy/het_30_tr_roseeta_data_dncon2_dist/train_list_0116221/test_list.txt'
# test_file = '/home/rajroy/deeplearning_codes/het_30_tr_roseeta_data_dncon2_divi_feature/training_list_het_121220/test_list.txt'
# test_file = '/home/rajroy/deeplearning_codes/het_30_tr_roseeta_data_dncon2_divi_feature/training_list_het_121220/test_list.txt'
test_file = '/home/rajroy/deeplearning_codes/het_30_tr_roseeta_data_dncon2_dist/train_list_0116221/test_list.txt'
# recall just extract how many
# cmap_dir = "/home/rajroy/test/cmap_divi_feature_583_l_200/"
# cmap_dir="/home/rajroy/q3_het30/some/"
cmap_dir="/home/rajroy/cmap_predict_400_hetero_f_582/"
test_file_name = file_reader(test_file)
val_array = []
all_threshold_values = []

report =1
# for n in range(5,100,5):
# THRESHOLD = n/100
SAMPLE_SIZE=0
for file in test_file_name:
    name = file.split('_')
    ##for distsss
    # predict_cmap = cmap_dir + file + '_.rr.npy.txt'
    # predict_cmap = cmap_dir + name[0]+'_'+name[1] + '.txt'
    ## for zdock
    # predict_cmap = cmap_dir + name[0]+'_'+name[1] + '.txt'

    ## for contacts
    predict_cmap = cmap_dir + name[0]+'_'+name[1] +'_'+str(name[2])+ '_.rr.npy.txt'
    if os.path.isfile(predict_cmap):
        SAMPLE_SIZE=SAMPLE_SIZE+1

        # /home/rajroy/predict_cmap_200_hetero/1H3OB_1H3OC_3305_.rr.npy.txt'
        ## for contacts
        real_cmap = cmap_dir + 'Y-' + file + '.txt.npy.txt'

        ##for distsss
        # real_cmap = cmap_dir  + file + '.txt.npy.txt'
        # /home/rajroy/predict_cmap_200_hetero/Y-1H3OB_1H3OC_3305.txt.npy.txt'

        if not os.path.exists(real_cmap):
            continue
        pred_arr = getY(predict_cmap)
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
        real_arr = getY(real_cmap)
        # total =math.sqrt(len_a*len_b)
        total = len_a + len_b
        # total = (len_a * len_b)**0.5
        # total =( len_a+ len_b)/2
        # rea_eval = np.where(real_arr == 1.0)
        # total =   len(rea_eval[0])
        new_pred_map = fix_pred_map(pred_arr, len_a, len_b)
        relax_0.append(calculateEvaluationStats(new_pred_map, real_arr, total,os.path.basename(file).split(".")[0]))

        real_arr_1 = make_relax(real_arr, 1)
        relax_1.append(calculateEvaluationStats(new_pred_map, real_arr_1, total,os.path.basename(file).split(".")[0]))

        real_arr_2 = make_relax(real_arr, 2)
        relax_2.append(calculateEvaluationStats(new_pred_map, real_arr_2, total,os.path.basename(file).split(".")[0]))

print(
    'RELAX ' + '\t\t\t' + 'TOP-5' + '\t\t\t' + 'TOP-10' + '\t\t\t' + 'TOP-20' + '\t\t\t' + 'TOP-30' + '\t\t\t' + 'TOP-50' + '\t\t\t' + 'L/30' + '\t\t\t' + 'L/20' + '\t\t\t' + 'L/10' +
    '\t\t\t' + 'L/5' + '\t\t\t' + 'L/' + '\t\t\t' + '2L/')
# relax_data_0=get_evaluation_result(relax_0,0,len(test_file_name))
# relax_data_1=get_evaluation_result(relax_1,1,len(test_file_name))
# relax_data_2=get_evaluation_result(relax_2,2,len(test_file_name))

relax_data_0=get_evaluation_result(relax_0,0,SAMPLE_SIZE)
relax_data_1=get_evaluation_result(relax_1,1,SAMPLE_SIZE)
relax_data_2=get_evaluation_result(relax_2,2,SAMPLE_SIZE)
print(SAMPLE_SIZE)


def report_individual_target(_data_array,_file_name):
    data_array=copy.deepcopy(_data_array)


    name_of_output_file =_file_name+ '.csv'

    with open(name_of_output_file, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(
            ['Name', 'TOP-T5', 'TOP-T10', 'TOP-T20', 'TOP-T30', 'TOP-T50', 'TOP-L30', 'TOP-L20', 'TOP-L10',   'TOP-L5', 'TOP-L', 'TOP-2L','TOTAL','RATIO'])
        for data in data_array:
            filewriter.writerow(data)
    print(output_dir + name_of_output_file)


if report == 1:

    output_dir = '/home/rajroy/'
    FILE_NAME = "het_30_578_feature_l_400_last_milstone_020121"
    report_individual_target(relax_data_0,output_dir+FILE_NAME+'_relax_0')
    report_individual_target(relax_data_1, output_dir + FILE_NAME + '_relax_1')
    report_individual_target(relax_data_2, output_dir + FILE_NAME + '_relax_2')