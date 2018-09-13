params.fasta = ""
params.references = ""
params.outdir = "results"
params.cpu = 4

fasta = Channel.from(file(params.fasta))
references = Channel.from(file(params.references))

process alignLocalRef {
  input:
  file fasta from fasta
  file references from references

  output:
  file "${fasta.simpleName}_ref.trim.aln" into aln

  publishDir params.outdir, mode: 'copy'
  tag "${fasta.simpleName}"
  cpus params.cpu

  script:
  """
  cat $fasta $references > ${fasta.simpleName}_ref.fasta
  mafft-einsi --adjustdirection --thread ${task.cpus} ${fasta.simpleName}_ref.fasta > ${fasta.simpleName}_ref.aln
  trimal -in ${fasta.simpleName}_ref.aln -out ${fasta.simpleName}_ref.trim.aln -gappyout
  """
}

process aln2Tree {
  input:
  file aln from aln

  output:
  file "${aln.simpleName}.treefile" into treefiles
  file "${aln.simpleName}.ufboot" into ufboots

  publishDir params.outdir, mode: 'copy'
  tag "${aln.simpleName}"
  cpus params.cpu

  script:
  """
  iqtree -s $aln -pre ${aln.simpleName} -bb 1000 -nt AUTO -ntmax ${task.cpus} -m TESTNEW -mset raxml -keep-ident -seed 12345 -wbtl
  """
}
