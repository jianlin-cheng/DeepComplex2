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


def calculateEvaluationStats(_pred_cmap, _true_cmap, L, _name):
    pred_cmap = copy.deepcopy(_pred_cmap)
    true_cmap = copy.deepcopy(_true_cmap)
    prec_T5, prec_T10, prec_T20, prec_T30, prec_T50, prec_L30, prec_L20, prec_L10, prec_L5, prec_L, prec_2L, con_num = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    max_Top = int((2 * L) + 0.5)
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
        if i == int((L / 2) + 0.5):
            prec_L2 = con_num * 100 / i
            if prec_L2 > 100: prec_L2 = 100
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((L) + 0.5):
            prec_L = con_num * 100 / i
            if prec_L > 100: prec_L = 100
            print("L=", L, "Val=", i, "Con_num=", con_num)
        if i == int((2 * L) + 0.5):
            prec_2L = con_num * 100 / i
            if prec_2L > 100: prec_2L = 100
            print("L=", L, "Val=", i, "Con_num=", con_num)

    return [prec_T5, prec_T10, prec_T20, prec_T30, prec_T50, prec_L30, prec_L20, prec_L10, prec_L5,prec_L2 ,prec_L, prec_2L,
            _name]


def get_evaluation_result(_arr, _relax, _SAMPLE_SIZE):
    sum_prec_T5, sum_prec_T10, sum_prec_T20, sum_prec_T30, sum_prec_T50, sum_prec_L30, sum_prec_L20, sum_prec_L10, sum_prec_L5,sum_prec_L2, sum_prec_L, sum_prec_2L = 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0
    # SAMPLE_SIZE = len(test_file_name)
    some_array = []
    # fasta_dict = loadFastaDictionary('/home/rajroy/Downloads/experiment_batch/fasta_dictionary.txt')
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
        sum_prec_L2 = sum_prec_L2 + values[9]
        sum_prec_L = sum_prec_L + values[10]
        sum_prec_2L = sum_prec_2L + values[11]

        # name_arrr = values[12].replace("__", "_").split('_')


        some_array.append(
            [values[12], str(values[0])[0:5], str(values[1])[0:5], str(values[2])[0:5], str(values[3])[0:5],
             str(values[4])[0:5],
             str(values[5])[0:5], str(values[6])[0:5], str(values[7])[0:5], str(values[8])[0:5], str(values[9])[0:5],
             str(values[10])[0:5]])

    print(str(_relax) + '\t\t\t' + str(sum_prec_T5 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_T10 / _SAMPLE_SIZE)[
                                                                                     0:5] + '\t\t\t' + str(
        sum_prec_T20 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_T30 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
        sum_prec_T50 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_L30 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
        sum_prec_L20 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_L10 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
        sum_prec_L5 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + '\t\t\t' + str(
        sum_prec_L2 / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(sum_prec_L / _SAMPLE_SIZE)[0:5] + '\t\t\t' + str(
        sum_prec_2L / _SAMPLE_SIZE)[0:5])

    return some_array



def rr2cmap(_rr,_cmap):
    rr_array =file_reader(_rr)[1:]
    for val in rr_array:
        values = val.split(" ")
        _cmap[int(values[0])-1][int(values[1])-1]=1
        _cmap[int(values[1]) - 1][int(values[0]) - 1] = 1
    return _cmap

def specific_filename_reader(_input_dir, _extension):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file:
                file_names.append(file.split(".")[0])
    return file_names


def fix_pred_map (_input):

    len = _input.shape[0]
    out = np.zeros((len, len))
    for i in range(len-1):
        for j in range(len - 1):
            temp = float((_input[i][j] + _input[j][i])) / 2
            out[i][j] = temp
            out[j][i] = temp
    return  out


relax_0 = []
relax_1 = []
relax_2 = []
# fasta_dict = loadFastaDictionary('/home/rajroy/Downloads/experiment_batch/fasta_dictionary.txt')
# test_file = '/home/rajroy/dncon2_inter_test_41.list'
# test_file = '/home/rajroy/predicted_contacts/list.txt'
# test_file = '//home/rajroy/test_set.txt'
# test_file = '/home/rajroy/test_list.txt'

# test_file_name = file_reader(test_file)
val_array = []
all_threshold_values = []
reject_list =["5YKX","4E1P","5YKZ","1IQ6"]
report = 1
# for n in range(5,100,5):
# THRESHOLD = n/100
SAMPLE_SIZE = 0

# predict_cmap_dir = "//home/rajroy/predicted_contacts//"
# real_cmap_dir = "//media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/hdd/DeepHomo/DeepHomo_testset/Altered_benchmark/labels_cmap/"
real_cmap_dir = "//home/rajroy/rr2cmap_elham/"
# real_cmap_dir = "//home/rajroy/Y-Labels/"
# real_rr =specific_filename_reader(real_cmap_dir,".cmap")
missing_list = []
# predict_cmap_dir=real_cmap_dir ="/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/hdd/DIMER_PAPER/deepHomo_casp/deepHomo_cmap/"
# predict_cmap_dir="/home/rajroy/experiment/cmap_symmetrical/"
# predict_cmap_dir="/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/hdd/DIMER_PAPER/homo_std_test_set/cmapfinal_deephomo_ep_82/"
predict_cmap_dir="/home/rajroy/predictions/"
real_rr =specific_filename_reader(predict_cmap_dir,".cmap")
for file in real_rr:
    print(file)
    if not file in reject_list :
        # name = file.split('_')
        predict_cmap =predict_cmap_dir+file +".cmap"
        # real_cmap = real_cmap_dir+file
        # real_cmap = real_cmap_dir+"Y-"+file.replace(".cmap",".txt")
        # real_cmap = real_cmap_dir+"Y-"+file+".txt"
        real_cmap = real_cmap_dir +file.replace("cmap_","")+".rr.cmap"

        # f = list(filter(lambda x: file in x, real_rr))
        # if len(f)==0:
        #     print("zero")
        #
        # # for values in f:
        #     # cmd_fasta = "cp " + values + " " + out_pdb + os.path.basename(values)
        # if len(f)>0:
        #     real_cmap = real_cmap_dir+"/"+f[0]+".txt"
        # else:
        #     continue
        # if not  os.path.exists(real_cmap):
        #     # real_cmap = real_cmap_dir + "/" + f[0] + ".rr"
        #     # real_cmap = real_cmap_dir + "/" + f[0]
        #     real_cmap = real_cmap_dir + "/" + "Y-"+file+".txt.npy.txt"

            # print(cmd_fasta)
            # os.system(cmd_fasta)
        # real_cmap = real_cmap_dir + "/" + "Y-" + file + ".txt"
        # predict_cmap = cmap_dir + name[0]+'_'+name[1] +'_'+str(name[2])+ '_.rr.npy.txt'
        if os.path.isfile(real_cmap) and os.path.isfile(predict_cmap):
            SAMPLE_SIZE=SAMPLE_SIZE+1
        else:
            missing_list.append(file)
            continue

        # if os.path.isfile(predict_cmap):
        #     # SAMPLE_SIZE = SAMPLE_SIZE + 1
        # else:
        #     missing_list.append(file)
        #     continue

        # real_cmap = cmap_dir + 'Y-' + file + '.txt.npy.txt'


        # if not os.path.exists(real_cmap):
        #     continue
        # pred_arr = getY(predict_cmap)
        # pred_arr = fix_pred_map(np.loadtxt(predict_cmap))
        pred_arr = np.loadtxt(predict_cmap)
        empty_cmap= np.zeros(pred_arr.shape)
        name = os.path.basename(predict_cmap).replace('.txt', '')
        # real_arr = rr2cmap(real_cmap,empty_cmap)
        real_arr = np.loadtxt(real_cmap)


        relax_0.append(calculateEvaluationStats(pred_arr, real_arr, real_arr.shape[0],os.path.basename(real_cmap).split(".")[0]))
        #
        real_arr_1 = make_relax(real_arr, 1)
        relax_1.append(calculateEvaluationStats(pred_arr, real_arr_1, real_arr.shape[0],os.path.basename(real_cmap).split(".")[0]))
        #
        real_arr_2 = make_relax(real_arr, 2)
        relax_2.append(calculateEvaluationStats(pred_arr, real_arr_2, real_arr.shape[0],os.path.basename(real_cmap).split(".")[0]))

        print(
            'RELAX ' + '\t\t\t' + 'TOP-5' + '\t\t\t' + 'TOP-10' + '\t\t\t' + 'TOP-20' + '\t\t\t' + 'TOP-30' + '\t\t\t' + 'TOP-50' + '\t\t\t' + 'L/30' + '\t\t\t' + 'L/20' + '\t\t\t' + 'L/10' +
            '\t\t\t' + 'L/5' + '\t\t\t' +'\t\t\t' + 'L/2' + '\t\t\t' + 'L/' + '\t\t\t' + '2L/')

    relax_data_0 = get_evaluation_result(relax_0, 0, SAMPLE_SIZE)
    relax_data_1 = get_evaluation_result(relax_1, 1, SAMPLE_SIZE)
    relax_data_2 = get_evaluation_result(relax_2, 2, SAMPLE_SIZE)
print(SAMPLE_SIZE)


def report_individual_target(_data_array, _file_name):
    data_array = copy.deepcopy(_data_array)

    name_of_output_file = _file_name + '.csv'

    with open(name_of_output_file, 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(
            ['Name', 'TOP-T5', 'TOP-T10', 'TOP-T20', 'TOP-T30', 'TOP-T50', 'TOP-L30', 'TOP-L20', 'TOP-L10', 'TOP-L5',
             'TOP-L', 'TOP-2L', 'TOTAL', 'RATIO'])
        for data in data_array:
            filewriter.writerow(data)
    print(output_dir + name_of_output_file)

print(missing_list)
if report ==1:
    output_dir = '/home/rajroy/'
    FILE_NAME = "homo_115_elham"
    report_individual_target(relax_data_0, output_dir + FILE_NAME + '_relax_0')
    report_individual_target(relax_data_1, output_dir + FILE_NAME + '_relax_1')
    report_individual_target(relax_data_2, output_dir + FILE_NAME + '_relax_2')
