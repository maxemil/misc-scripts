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
overlap = Channel.from("55 60 65 70 75 80 85 90".tokenize())
ident = Channel.from("05 10 15 20 25 30 35 40".tokenize())

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

  script:
  """
  #!/usr/bin/env anapy3

  import pandas as pd
  from math import isclose

  prefix = "${fnodes.simpleName}"
  try:
      print("read {}.fnodes".format(prefix))
      fnodes = pd.read_csv("$fnodes", header=None, sep='\\t|&', engine='python')
      print("read {}.cluster".format(prefix))
      cluster = pd.read_csv("$cluster", header=None, names=["clst"])
      print("create members column for {}".format(prefix))
      cluster['members'] = cluster.apply(lambda row: (fnodes[0] == row['clst']).sum(), axis=1)
      print("create spec column for {}".format(prefix))
      cluster['spec'] = cluster.apply(lambda row: len(set(fnodes[1][fnodes[0] == row['clst']])), axis=1)
      print("create panorthologs file for {}".format(prefix))
      panorthologs = cluster[cluster.apply(lambda row: isclose(row['members'], row['spec'], abs_tol=5) and row['spec'] > (0.9 * $num_taxa), axis=1)]
      panorthologs.to_csv("{}.panorthologs".format(prefix), sep='\\t', index=False)
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


  panorthologs = glob.glob("*.panorthologs")
  pan_num = pd.DataFrame(columns=["05", "10", "15", "20", "25", "30", "35", "40"], index=["55", "60", "65", "70", "75", "80", "85", "90"])

  for p in panorthologs:
      ident = p.replace('.panorthologs', '').split('_')[1].replace('i', '')
      overlap = p.replace('.panorthologs', '').split('_')[2].replace('o', '')
      pan_list = pd.read_csv(p, sep='\\t', header=1)
      pan_num[ident][overlap] = pan_list.shape[0]

  pan_num = pan_num.fillna(0)
  ax = sns.heatmap(pan_num, linewidth=0.5, annot=True, fmt="d", cbar=False)
  ax.set(xlabel="identity", ylabel="overlap")
  plt.savefig('panorthologs_numbers.pdf')

  max_value = pan_num.max(1).max()
  for k,v in dict(pan_num.idxmax()).items():
      if pan_num[k][v] == max_value:
          print("silix_i{}_o{}".format(k,v), end="")
  """
}
silix_optimal = silix_get_optimal.spread(x).filter{it[0].simpleName == it[3]}


process HiFixClustering {
  input:
  set file(fnodes), file("allTaxa.fasta"), file(net), val(prefix) from silix_optimal

  output:
  file "${fnodes.simpleName}.hifix.fnodes" into hifix_fnodes
  file "${fnodes.simpleName}.hifix.fsize" into hifix_fsize
  file "*.fasta" into cluster_fasta

  container "HiFiX.img"
  publishDir "$params.outdir/hifix", mode: 'copy'
  tag "$prefix"

  script:
  """
  hifix -t 20 -n 51 allTaxa.fasta $net $fnodes > ${fnodes.simpleName}.hifix.fnodes
  silix-split -n 1 allTaxa.fasta ${fnodes.simpleName}.hifix.fnodes -p $prefix
  silix-fsize ${fnodes.simpleName}.hifix.fnodes > ${fnodes.simpleName}.hifix.fsize
  rename "allTaxa_silix_i30_o55" "" *.fasta
  """
}
