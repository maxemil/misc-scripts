params.cluster = ""
params.output_alignments = "alignments"
params.output_faa = "cluster_faa"

cluster = Channel.fromPath(params.cluster)

process sortClustersBySize {
  input:
  file faa from cluster

  output:
  file "small_cluster/${faa.simpleName}.small.faa" into small_faa optional true
  file "single_cluster/${faa.simpleName}.single.faa" into single_faa optional true
  file "${faa.simpleName}.faa" into normal_faa optional true

  publishDir "${params.output_faa}", mode: 'copy'
  tag {"${faa.simpleName}"}

  script:
  """
  if [ "\$(grep -c '^>' $faa)" -lt 4 ] && [ "\$(grep -c '^>' $faa)" -gt 1 ]
  then
    mkdir small_cluster
    mv $faa small_cluster/${faa.simpleName}.small.faa
  elif [ "\$(grep -c '^>' $faa)" -eq 1 ]
  then
    mkdir single_cluster
    mv $faa single_cluster/${faa.simpleName}.single.faa
  else
    cp -L $faa ${faa.simpleName}.faa
  fi
  """
}

process alignCluster {
  input:
  file faa from normal_faa

  output:
  file "${faa.baseName}.aln" into trimmed_alignments

  publishDir "${params.output_alignments}", mode: 'copy'
  tag {"${faa.simpleName}"}
  cpus 1

  script:
  """
  sed -i 's/*//g' $faa
  prequal $faa
  mafft-einsi --anysymbol --thread ${task.cpus} ${faa}.filtered > "${faa.baseName}.aln"
  #trimal -in "${faa.baseName}.aln" -out "${faa.baseName}.trim.aln" -gappyout -fasta
  """
}
