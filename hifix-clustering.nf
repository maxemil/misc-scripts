params.faa = ""
params.outdir = "results"
params.cpu = 4

Channel.from(file(params.faa)).into{count_faa; merge_faa}
merge_faa.collectFile(name: 'allTaxa.faa', storeDir: params.outdir).into{subject_faa; query_faa; silix_faa}
num_taxa = count_faa.count()

process diamondMakeDBFromFaas {
  input:
  file faa from subject_faa

  output:
  file "${faa.simpleName}.dmnd" into database

  tag {"${faa.simpleName}"}

  script:
  """
  diamond makedb --in $faa --db ${faa.simpleName}.dmnd
  """
}


process DiamondBlastAllvsAll {
  input:
  file faa from query_faa
  file database from database

  output:
  file "${faa.simpleName}.blastp" into blast_results
  file "${faa.simpleName}.unaligned.blastp" into blast_unaligned optional true

  publishDir params.outdir, mode: 'copy'
  tag "${faa.simpleName}"
  cpus params.cpu

  script:
  """
  diamond blastp --query $faa \
               --db $database \
               --out ${faa.simpleName}.blastp \
               --threads ${task.cpus} \
               --outfmt 6 \
               --max-target-seqs 1000 \
               --more-sensitive \
               --evalue 0.0001 \
  """
}
overlap_array = "55 60 65 70 75 80 85 90".tokenize()
// overlap_array = "55 60 65 70 75".tokenize()
overlap = Channel.from(overlap_array)
ident_array = "05 10 15 20 25 30 35 40".tokenize()
// ident_array = "20 25 30 35".tokenize()
ident = Channel.from(ident_array)

silix_params = overlap.spread(ident)

process SilixClustering {
  input:
  file blastp from blast_results.first()
  file faa from silix_faa.first()
  set val(o), val(i) from silix_params

  output:
  set file("silix_i${i}_o${o}.fnodes"), file("silix_i${i}_o${o}.cluster") into silix_results
  set file("silix_i${i}_o${o}.fnodes"), file("$faa"), file("${faa.simpleName}.net") into silix_get_optimal

  publishDir "$params.outdir/silix", mode: 'copy'
  stageInMode "copy"
  container "HiFiX.img"
  tag "silix i=0.$i o=0.$o"

  script:
  """
  realpath $blastp > blastfilepath
  silix $faa blastfilepath --prefix silix_i${i}_o${o} --ident 0.$i --overlap 0.$o --net > silix_i${i}_o${o}.fnodes
  cut -f1 silix_i${i}_o${o}.fnodes | sort | uniq  > silix_i${i}_o${o}.cluster
  """
}


process CountSilixPanorthologs {
  input:
  set file(fnodes), file(cluster) from silix_results
  val num_taxa from num_taxa
  output:
  file "${fnodes.simpleName}.panorthologs" into panorthologs

  publishDir "$params.outdir/silix", mode: 'copy'
  tag "${fnodes.simpleName}"
  cpus params.cpu

  script:
  """
  #!/usr/bin/env anapy3

  import pandas as pd
  from math import isclose
  import numpy as np
  from multiprocessing import Pool
  from functools import partial

  def count_panorthologs(row, fnodes):
      spec = len(set(fnodes[1][fnodes[0] == row['clst']]))
      members = (fnodes[0] == row['clst']).sum()
      if isclose(members, spec, abs_tol=5) and spec > (0.9 * $num_taxa): # count universal single copy gene families
      # if spec > (0.9 * $num_taxa): # count universal gene families
          return True
      else:
          return False

  def count_panorthologs_df(df, fnodes):
      counts = df.apply(partial(count_panorthologs, fnodes=fnodes), axis=1)
      print(counts)
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
  """
}


process PlotSilixPanorthologs {
  input:
  file "*" from panorthologs.collect()

  output:
  file 'panorthologs_numbers.pdf' into panortholog_plot
  stdout x

  publishDir "$params.outdir/silix", mode: 'copy'
  tag "all"

  script:
  """
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
  """
}
silix_optimal = silix_get_optimal.spread(x).filter{it[0].simpleName == it[3]}


process HiFixClustering {
  input:
  set file(fnodes), file("allTaxa.fasta"), file(net), val(prefix) from silix_optimal
  val num_taxa from num_taxa

  output:
  file "${fnodes.simpleName}.hifix.fnodes" into hifix_fnodes
  file "${fnodes.simpleName}.hifix.fsize" into hifix_fsize
  file "*.fasta" into cluster_fasta

  container "HiFiX.img"
  publishDir "$params.outdir/hifix", mode: 'copy'
  tag "$prefix"
  cpus params.cpu

  script:
  """
  hifix -t ${task.cpus} -n $num_taxa --force allTaxa.fasta $net $fnodes > ${fnodes.simpleName}.hifix.fnodes
  silix-split -n 1 allTaxa.fasta ${fnodes.simpleName}.hifix.fnodes -p $prefix
  silix-fsize ${fnodes.simpleName}.hifix.fnodes > ${fnodes.simpleName}.hifix.fsize
  find . -name '*.fasta' -exec rename "allTaxa_$prefix" "" {} +
  """
}
