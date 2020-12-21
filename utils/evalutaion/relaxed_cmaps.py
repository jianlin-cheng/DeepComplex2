#takes in a cmap and gets the relax maps
import copy
import os

import numpy as np

input_cmap = '/home/rajroy/1A6ZA_1A6ZB.txt'


def make_relax(_Y, _relax):
    comb = []
    if _relax == 1:
        comb = [-1, 0, 1]
    if _relax == 2:
        comb = [-2, -1, 0, 1, 2]

    original_arr = copy.deepcopy(_Y)

    result = np.where(_Y == 1.0)
    number_contacts = len(result[0])
    # print('before relax ' + str(number_contacts))
    counter_tracker = 0
    for val in range(0, number_contacts):
        for x_increment in comb:
            for y_increment in comb:
                index_X = result[0][val]
                index_Y = result[1][val]
                # be concerned wheter the alter indexes will be greater or less than the legit value and hence checking
                if index_X + x_increment >= 0 and index_Y + y_increment >= 0:
                    if index_X + x_increment < len(original_arr) and index_Y + y_increment < len(original_arr):
                        # print(index_X + x_increment)
                        # print(index_Y + y_increment)
                        counter_tracker = counter_tracker + 1
                        original_arr[index_X + x_increment][index_Y + y_increment] = 1.0
                        # print(original_arr[index_X + x_increment][index_Y + y_increment])
                        # print(original_arr[index_X][index_Y])
    result = np.where(original_arr == 1.0)
    number_contacts = len(result[0])
    # print('afterrelax ' + str(number_contacts))
    # print(counter_tracker)

    return original_arr



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


# relax_0_array= getY(input_cmap)
# relax_1_array = make_relax(relax_0_array, 1)
# relax_2_array = make_relax(relax_0_array, 2)
# # print(y)
