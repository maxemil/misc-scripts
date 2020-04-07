params.cluster = ""
params.output_alignments = "alignments"
params.output_faa = "cluster_faa"

cluster = Channel.fromPath(params.cluster)

process sortClustersBySize {
  input:
  file faa from cluster

  output:
  file "small_cluster/${faa.simpleName}.faa" into small_faa optional true
  file "${faa.simpleName}.faa" into normal_faa optional true

  publishDir "${params.output_faa}", mode: 'copy', pattern: 'small_cluster/*'
  tag {"${faa.simpleName}"}

  script:
  """
  if [ "\$(grep -c '^>' $faa)" -lt 4 ]
  then
    mkdir small_cluster
    mv $faa small_cluster/${faa.simpleName}.faa
  else
    cp -L $faa ${faa.simpleName}.faa
  fi
  """
}

normal_faa.into{normal_faa_alignment; normal_faa_filering}


process alignCluster {
  input:
  file faa from normal_faa_alignment

  output:
  file "${faa.simpleName}.divvy.aln" into divvyied_alignments
  file "${faa.simpleName}.aln" into untreated_alignments

  // publishDir "${params.output_alignments}", mode: 'copy'
  tag {"${faa.simpleName}"}
  cpus 30

  script:
  """
  sed -i 's/*//g' $faa
  prequal $faa
  mafft-einsi --anysymbol --thread ${task.cpus} ${faa}.filtered > "${faa.simpleName}.aln"
  divvier -divvy -divvygap "${faa.simpleName}.aln"
  mv ${faa.simpleName}.aln.divvy.fas ${faa.simpleName}.divvy.aln
  """
}

filter_alignments = normal_faa_filering.combine(divvyied_alignments)

process filterAlignmentSize{
  input:
  set file(faa), file(divvy) from filter_alignments

  output:
  file "$divvy" into filtered_divvyied_alignments optional true
  file "$faa" into filtered_faa optional true
  file "small_cluster/*.faa" into filtered_small_cluster optional true

  publishDir "${params.output_alignments}", mode: 'copy', pattern: '*.aln'
  publishDir "${params.output_faa}", mode: 'copy', pattern: 'small_cluster/*'
  publishDir "${params.output_faa}", mode: 'copy', pattern: '*.faa'
  stageInMode 'copy'
  tag {"${faa.simpleName}"}

  when:
  "${faa.simpleName}" == "${divvy.simpleName}"

  script:
  """
  #! /usr/bin/env python3
  from Bio import SeqIO
  import os

  os.makedirs("small_cluster")
  missing_data = ['-', 'X']

  divvy = {rec.id:rec for rec in SeqIO.parse("$divvy", 'fasta')}
  faa = {rec.id:rec for rec in SeqIO.parse("$faa", 'fasta')}
  all_missing_data = []
  for recid, rec in divvy.items():
        if all([c in missing_data for c in rec.seq]):
            all_missing_data.append(recid)
  [divvy.pop(recid) for recid in all_missing_data]
  small = "" if len(divvy.keys()) > 3 else 'small_cluster/'
  with open("{}$faa".format(small), 'w') as faa_file, open("{}$divvy".format(small), 'w') as divvy_file:
      for recid, rec in faa.items():
          if recid in divvy.keys():
              SeqIO.write(divvy[recid], divvy_file, 'fasta')
              SeqIO.write(rec, faa_file, 'fasta')
          else:
              with open("small_cluster/{}.faa".format(recid.replace('.', '_')), 'w') as single:
                  SeqIO.write(rec, single, 'fasta')
  if small:
    os.remove("$divvy")
    os.remove("$faa")
  """

}
