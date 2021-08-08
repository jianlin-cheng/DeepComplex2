import copy

import numpy as np

def bin_to_label(_input):
    return  np.sum(_input[1:7])



def dist_to_label(_input,_len_a,_len_b):
    value = copy.deepcopy(_input)
    total = _len_a + _len_b
    label_dist = np.zeros((total,total))
    for j_counter in range(0, total):
        for i_counter in range(0, total):

            label_dist[j_counter][i_counter] = bin_to_label(value[j_counter][i_counter] )

    print(label_dist)
    #
    # for j_counter in range(0, _len_a):
    #     for i_counter in range(0, _len_a):
    #         output_rr[j_counter][i_counter] = 0
    #         output_rr[i_counter][j_counter] = 0
    #
    # for j_counter in range(_len_a, total):
    #     for i_counter in range(_len_a, total):
    #         output_rr[j_counter][i_counter] = 0
    #         output_rr[i_counter][j_counter] = 0
    #
    #
    # for j_counter in range(0, total):
    #     for i_counter in range(0, total):
    #         value[j_counter][i_counter] = output_rr[i_counter][j_counter] + output_rr[j_counter][i_counter]
    #         value[i_counter][j_counter] = 0


content = np.load('/media/rajroy/fbc3794d-a380-4e0f-a00a-4db5aad57e75/rajroy/het_30_bin_200/1A15A_1A15B.npz')
y_labels =content.f.arr_0.squeeze()
print(y_labels)
dist_to_label(y_labels,67,57)