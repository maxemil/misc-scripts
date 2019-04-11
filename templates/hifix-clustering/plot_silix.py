#!/usr/bin/env anapy3

import pandas as pd
import glob
import seaborn as sns
import matplotlib.pylab as plt
from itertools import product
import numpy as np

ident_list = "$ident_array".replace('[', '').replace(']', '').replace(' ', '').split(',')
overlap_list = "$overlap_array".replace('[', '').replace(']', '').replace(' ', '').split(',')
panorthologs = glob.glob("*.panorthologs")
pan_num = pd.DataFrame(columns=ident_list, index=overlap_list)

for p in panorthologs:
  ident = p.replace('.panorthologs', '').split('_')[1].replace('i', '')
  overlap = p.replace('.panorthologs', '').split('_')[2].replace('o', '')
  pan_list = pd.read_csv(p, sep='\\t', header=1)
  pan_num[ident][overlap] = pan_list.shape[0]
pan_num = pan_num.fillna(0)

max_value = pan_num.max(1).max()
maxima = []
for i,o in product(ident_list, overlap_list):
  if pan_num[i][o] == max_value:
      maxima.append((i,o))

maximum = (0,0,0)
for m in maxima:
  col_rowsum = sum(pan_num[m[0]]) + sum(pan_num.loc[m[1]])
  if col_rowsum > maximum[2]:
      maximum = (m[0], m[1], col_rowsum)

print("silix_i{}_o{}".format(maximum[0],maximum[1]), end="")

mask_array = np.zeros(pan_num.shape, dtype=bool)
mask_array[pan_num.index.get_loc(maximum[1])][pan_num.columns.get_loc(maximum[0])] = True
mask_array = mask_array == False
ax = sns.heatmap(pan_num, linewidth=0.5, annot=True, fmt="d", cbar=False)
ax = sns.heatmap(pan_num, linewidth=0.5, annot=True, fmt="d", cbar=False, annot_kws={'weight': 'bold', 'color': 'red'}, mask=mask_array)
ax.set(xlabel="identity", ylabel="overlap")
plt.savefig('panorthologs_numbers.pdf')
