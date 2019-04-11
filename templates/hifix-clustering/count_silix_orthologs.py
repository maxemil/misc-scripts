#!/usr/bin/env anapy3

import pandas as pd
from math import isclose
import numpy as np
from multiprocessing import Pool
from functools import partial

def count_panorthologs(row, fnodes):
  spec = len(set(fnodes[1][fnodes[0] == row['clst']]))
  members = (fnodes[0] == row['clst']).sum()
  # if spec > (0.9 * $num_taxa): # count universal gene families
  if isclose(members, spec, abs_tol=5) and spec > (0.9 * $num_taxa): # count universal single copy gene families
      return True
  else:
      return False

def count_panorthologs_df(df, fnodes):
  counts = df.apply(partial(count_panorthologs, fnodes=fnodes), axis=1)
  return counts

def parallelize(cluster, fnodes, func, processors):
  with Pool(processors) as pool:
      data_split = np.array_split(cluster, processors)
      result_split = pool.map(partial(func, fnodes=fnodes), data_split)
      return pd.concat(result_split)

def create_panorthologs_files(prefix):
  print("read {}.fnodes".format(prefix))
  fnodes = pd.read_csv("$fnodes", header=None, sep='\\t|&', engine='python')
  print("read {}.cluster".format(prefix))
  cluster = pd.read_csv("$cluster", header=None, names=["clst"])
  print("get panorthologs counts for {}".format(prefix))
  panorthologs_count = parallelize(cluster, fnodes, count_panorthologs_df, ${task.cpus})
  print("extract panortholog clusters for {}".format(prefix))
  panorthologs = cluster[panorthologs_count]
  print("print panorthologs file for {}".format(prefix))
  panorthologs.to_csv("{}.panorthologs".format(prefix), sep='\\t', index=False)

if __name__ == "__main__":
  prefix = "${fnodes.simpleName}"
  try:
      create_panorthologs_files(prefix)
  except:
      print("error when processing {}".format(prefix))
