params.faa = ""
params.outdir = "results"
params.cpu = 4

Channel.fromPath(params.faa).into{count_faa; merge_faa}
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
  template 'hifix-clustering/count_silix_orthologs.py'
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
  template 'hifix-clustering/plot_silix.py'
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
  hifix -t ${task.cpus} -n 2 --force allTaxa.fasta $net $fnodes > ${fnodes.simpleName}.hifix.fnodes
  silix-split -n 1 allTaxa.fasta ${fnodes.simpleName}.hifix.fnodes -p $prefix
  silix-fsize ${fnodes.simpleName}.hifix.fnodes > ${fnodes.simpleName}.hifix.fsize
  find . -name '*.fasta' -exec rename "allTaxa_$prefix" "" {} +
  """
}
