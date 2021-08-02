import os

from numpy import genfromtxt

# result_csv = genfromtxt('///home/rajroy/regressional_plot/Aveerage_cmap_ep_82_testset_true_relax_0.csv', delimiter=',', dtype='str')
result_csv = genfromtxt('//home/rajroy/regressional_plot/ep_82_testset_predicted_relax_0.csv', delimiter=',', dtype='str')
# result_csv = genfromtxt('//home/rajroy/regressional_plot/intra_test_set_precision.csv', delimiter=',', dtype='str')
print(result_csv)


def specific_filename_reader(_input_dir, _extension):
    file_names = []
    for root, directories, files in os.walk(_input_dir):
        for file in files:
            if _extension in file:
                file_names.append(file.split(".")[0])
    return file_names


l5_array = {}
for val in result_csv[1:len(result_csv)]:
# for val in result_csv[1:50]:
#inter
    name = val[0].replace('[', "")
    name = name.replace(']', "")
    name = name.replace("'", "")
    # print(name)
    l5_array[name] = float(val[9])
    #intra
    # name=val[0].replace(".rr_0.5 (precision)","")
    # l5_array[name] = float(val[6])


import matplotlib.pyplot as plt;

plt.rcdefaults()
import numpy as np

inp_file = "//home/rajroy/Y-Labels/"
contact_info = specific_filename_reader(inp_file, '.txt')

O_5 = []
O_5_1 = []
O_1_15 = []
O_15_20 = []
O_20_25 = []
O_25_30 = []
O_30_35 = []
O_35_40 = []
# O_40_45 = []
# O_45_50 = []
counter = 0
for values in l5_array:
    counter =counter+1
    print(counter)
    # file_name = inp_file + "Y-" + values + ".txt"
    file_name = inp_file + values + ".txt"
    if not os.path.exists(file_name):
        continue
    val = np.loadtxt(file_name)
    length = val.shape[0]
    contact = np.count_nonzero(val)
    dif = float(contact / length)
    if dif <= 0.5:
        O_5.append(l5_array.get(values))
    elif dif <= 1.0:
        O_5_1.append(l5_array.get(values))
    elif dif <= 1.5:
        O_1_15.append(l5_array.get(values))
    elif dif <= 2.0:
        O_15_20.append(l5_array.get(values))
    elif dif <= 2.5:
        O_20_25.append(l5_array.get(values))
    elif dif <= 3.0:
        O_25_30.append(l5_array.get(values))
    elif dif <= 3.5:
        O_30_35.append(l5_array.get(values))
    elif dif <= 4.0:
        O_35_40.append(l5_array.get(values))
    # elif dif <= 4.5:
    #     O_40_45.append(l5_array.get(values))
    # elif dif <= 5.0:
    #     O_45_50.append(l5_array.get(values))

    # print(values)

    # append all precision based on density

# here give the average
def Average(lst):
    if       sum(lst)    ==0:
        return 0
    if len(lst) == 0:
        return 0
    return round(sum(lst) / len(lst), 5)

objects = ('0-0.5', '0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5', '2.5-3.0','3.0-3.5','3.5-4.0' )
# objects = ('0-0.5', '0.5-1.0', '1.0-1.5', '1.5-2.0', '2.0-2.5', '2.5-3.0','3.0-3.5','3.5-4.0','4.0-4.5','4.5-5.0')
y_pos = np.arange(len(objects))
# performance = [Average( O_5), Average(O_5_1 ),Average(O_1_15 ),Average(O_15_20 ),Average(O_20_25 ),Average(O_25_30 ),Average(O_30_35),Average( O_35_40 ),Average(O_40_45 ),Average(O_45_50 )]
performance = [Average( O_5), Average(O_5_1 ),Average(O_1_15 ),Average(O_15_20 ),Average(O_20_25 ),Average(O_25_30 ),Average(O_30_35),Average( O_35_40 )]




plt.bar(y_pos, performance, align='center', alpha=0.5)
plt.xticks(y_pos, objects)
# plt.ylabel('Average IntraChain Contact Precision for 2L in %')
# plt.title('Range of True Contact Densities')
plt.ylabel('Average Interchain Contact Precision for L/5 ( % )'  )
# plt.ylabel('Average Intrachain Contact Precision for 2L ( % )' )
# plt.title('Range of Interchain Contact Densities')
plt.xlabel('Range of Interchain Contact Densities')
plt.show()
