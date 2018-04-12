#!/usr/bin/env nextflow

params.faas = ""
params.outdir = ""
params.cpus = 2

faas = Channel.fromPath(params.faas)

process alignTrimRefQueries {
  input:
  file faa from faas

  output:
  file "${faa.baseName}.trim.aln" into trimmed_alignments
  file "${faa.baseName}.aln" into alignments

  publishDir "${params.outdir}/alignments", mode: 'copy'
  tag {"${faa.simpleName}"}

  script:
  """
  mafft-einsi --thread ${task.cpus} $faa > "${faa.simpleName}.aln"
  trimal -in "${faa.simpleName}.aln" -out "${faa.simpleName}.trim.aln" -gappyout -fasta
  """
}

process buildTrees {
  input:
  file aln from trimmed_alignments

  output:
  file "${aln.simpleName}.treefile" into treefiles
  file "${aln.simpleName}.ufboot" into ufboots

  container 'iqtree.simg'
  cpus "${params.cpus}"
  publishDir "${params.outdir}/ufboots", mode: 'copy'
  tag {"${aln.simpleName}"}

  script:
  """
  iqtree -s $aln \
             -m TEST \
             -mset LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60,EX2,EX3,EHO,UL2,UL3,EX_EHO,LG4M,LG4X \
             -wbt \
             -bb 1000 \
             -bnni \
             -nt AUTO \
             -ntmax ${task.cpus} \
             -pre ${aln.simpleName}
  """
}


process computeSplits {
  input:
  file ufboot from ufboots

  output:
  set val("${ufboot.simpleName}"), file("${ufboot.simpleName}.splits") into splits

  container 'discordance.simg'
  publishDir "${params.outdir}/splits", mode: 'copy'
  tag {"${ufboot.simpleName}"}

  script:
  """
  tre_make_splits.pl --file ${ufboot.simpleName}.ufboot --outfile ${ufboot.simpleName}.splits
  sed -i -e 's/\\.\\.[a-zA-Z0-9\\.\\_]*//g' ${ufboot.simpleName}.splits
  """
}

splits.into{splits1; splits2}

split_pairs = splits1.combine(splits2)

process discordancePair {
  input:
  set val(split1), file("splitfile1"), val(split2), file("splitfile2") from split_pairs

  output:
  file "${split1}.${split2}.out" into discordance_scores

  container 'discordance.simg'
  publishDir "${params.outdir}/discordance/pairwise", mode: 'copy'
  tag {"${split1}_${split2}"}

  when:
  "${split1}" != "${split2}"

  script:
  """
  tre_discordance_two.pl --file1 splitfile1 --file2 splitfile2 --outfile ${split1}.${split2}.out
  """
}

discordance_files = discordance_scores.collectFile(storeDir: "${params.outdir}/discordance") { item -> [ "${item.simpleName}.discordance", item ]}.collect()


process plotDiscordance {
  input:
  file all_scores from discordance_files

  output:
  file "discordance_score.pdf" into discordance_plot
  publishDir "${params.outdir}/discordance", mode: 'copy'

  script:
  """
  #!/usr/bin/env anapy3

  import pandas as pd
  import glob
  import seaborn as sns
  import matplotlib.pyplot as plt

  score_dict = {}
  scores = glob.glob("*.discordance")
  for handle in scores:
      fam = handle.replace('.discordance', '')
      score_dict[fam] = sum([float(score) for score in open(handle)])

  df = pd.DataFrame.from_dict(score_dict, orient='index')
  df = df.sort_values(by=0)
  df.index = range(1, df.size + 1)
  df.columns = ['summed score']
  ax = df.plot()
  ax.set_ylabel('Sum of pairwise discordance scores')
  ax.set_xlabel('OGs ranked by increasing discordance')

  plt.axvline(x=.7*df.size)
  plt.axvline(x=.9*df.size, ls='dashed')

  plt.savefig("discordance_score.pdf")
  """
}
