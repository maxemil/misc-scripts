#!/usr/bin/env nextflow

params.fasta = ""
params.outdir = ""
params.server_adress = ""

fasta_input = Channel.fromPath(params.fasta)

process runEMapper {
  input:
  file fasta from fasta_input

  output:
  file "${fasta.simpleName}.emapper.hmm_hits" into hmm_hits
  file "${fasta.simpleName}.emapper.seed_orthologs" into seed_orthologs
  file "${fasta.simpleName}.emapper.annotations" into annotations

  publishDir "${params.outdir}", mode: 'copy'
  tag {"${fasta.baseName}"}
  stageInMode 'copy'

  script:
  """
  python2 /local/two/Software/eggnog-mapper/emapper.py -d ${params.server_adress} -o ${fasta.simpleName} --output_dir . -i $fasta
  """
}
