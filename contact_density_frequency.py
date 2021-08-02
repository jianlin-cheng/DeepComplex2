import os

import numpy as np
from matplotlib import pyplot as plt
from numpy import genfromtxt


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
performance = [667,747,191,53,23,14,3,0]




plt.bar(y_pos, performance, align='center', alpha=0.5)
plt.xticks(y_pos, objects)
# plt.ylabel('Average IntraChain Contact Precision for 2L in %')
# plt.title('Range of True Contact Densities')
plt.ylabel('Frequency'  )
# plt.ylabel('Average Intrachain Contact Precision for 2L ( % )' )
# plt.title('Range of Interchain Contact Densities')
plt.xlabel('Range of Interchain Contact Densities')
plt.show()
