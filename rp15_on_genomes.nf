
params.genomes = ""
params.references = ""
params.faas = ""
params.phylo_method = "fasttree"
params.outdir = "outdir"
params.numprots = 5
params.cpus = 10
params.prune_percentage = "0.0"

params.get_hits = false

if(params.get_hits){
  Channel.from(params.get_hits).into{get_hits_database; get_hits_references}
}else{
  Channel.empty().into{get_hits_database; get_hits_references}
}

if(params.faas != ""){
  additional_faas = Channel.fromPath(params.faas)
}else{
  additional_faas = Channel.empty()
}
genomes = Channel.fromPath(params.genomes)
Channel.fromPath(params.references).into{references; references_add}

process predictGenes {
  input:
  file genome from genomes

  output:
  file "PROKKA/${genome.baseName}.faa" into faas
  file "PROKKA/${genome.baseName}.*" into prokka_output

  publishDir "${params.outdir}", mode: 'copy', overwrite: 'true'
  tag {"${genome.simpleName}"}

  script:
  """
  /local/two/Software/prokka-partial/bin/prokka --outdir PROKKA \
            --partialgenes --prefix ${genome.baseName} \
            --locustag ${genome.baseName} $genome
  """
}

faa_combined = faas.concat(additional_faas).collectFile(name: 'all.faa', newLine: true)

process databaseFromFaas {
  input:
  file all from faa_combined
  val get_hits from get_hits_database.first()

  output:
  file "${all.baseName}.db*" into database

  stageInMode 'copy'
  tag {"${all.baseName}"}

  script:
  """
  mv $all ${all.baseName}.db
  makeblastdb -in ${all.baseName}.db -dbtype prot -parse_seqids
  """
}

process alignTrimReferences {
  input:
  file ref from references
  val get_hits from get_hits_references.first()

  output:
  set file(ref), file("${ref.baseName}.trim.aln") into trimmed_references

  tag {"${ref.simpleName}"}

  script:
  """
  mafft-linsi --thread ${task.cpus} $ref > "${ref.baseName}.aln"
  trimal -in "${ref.baseName}.aln" -out "${ref.baseName}.trim.aln" -gappyout -fasta
  """
}

database.into{psiDB; retrieveDB}

process psiBlast {
  input:
  set file(ref), file(trim_ref) from trimmed_references
  file db from psiDB.first()

  output:
  set file("${ref.baseName}_vs_all.psiblast"), file(ref) into psiblast_results

  publishDir "${params.outdir}/psiblast"
  tag {"${ref.simpleName}"}

  script:
  """
  psiblast -in_msa $trim_ref \
    -db all.db \
    -num_alignments 1400 \
    -evalue 1e-6 \
    -outfmt '6 std qcovhsp stitle' \
    -out "${ref.baseName}_vs_all.psiblast" \
    -num_threads ${task.cpus}
  """
}

process addHits {
  input:
  file db from retrieveDB.first()
  set  file(psiblast), file(ref) from psiblast_results
  // file ref from references_add

  output:
  file "${ref.baseName}_added.faa" into added_faas

  tag {"${ref.simpleName}"}

  script:
  """
  cp -L $ref ${ref.baseName}_added.faa
  blastdbcmd \
    -db all.db -dbtype prot \
    -entry_batch <( awk '\$13 >= 70 && \$3 >= 25 {print \$0}' $psiblast | cut -f 2 | sed -r "s/lcl\\|//" | sort -u ) \
    -outfmt %f >> ${ref.baseName}_added.faa
    """
}

process alignTrimRefQueries {
  input:
  file faa from added_faas

  output:
  file "${faa.baseName}.trim.aln" into trimmed_alignments

  publishDir "${params.outdir}", mode: 'copy'
  tag {"${faa.simpleName}"}

  script:
  """
  mafft-linsi --thread ${task.cpus} $faa > "${faa.baseName}.aln"
  trimal -in "${faa.baseName}.aln" -out "${faa.baseName}.trim.aln" -gappyout -fasta
  """
}

all_trimmed_alignments = trimmed_alignments.collect()

process concatenateAlignments {
  input:
  file "*" from all_trimmed_alignments

  output:
  file "concat.fasta" into concatenated_alignment

  publishDir "${params.outdir}", mode: 'copy'
  stageInMode 'copy'

  script:
  """
  sed -i -E '/^>lcl/ s/_[0-9]*\$//g' *.aln
  sed -i -E '/^>lcl/! s/_gi\\|.*\$//g' *.aln
  sed -i -E '/^>lcl/! s/_[0-9]*_[3|5]0S_.*\$//g' *.aln
  sed -i -E '/^>lcl/! s/_intID[0-9]*\$//g' *.aln
  anapy3 $workflow.projectDir/scripts/concatenate.py *.aln -t $params.numprots
  """
}

process pruneAlignment {
  input:
  file aln from concatenated_alignment

  output:
  file "${aln.baseName}.pruned.aln" into pruned_alignment

  publishDir "${params.outdir}", mode: 'copy'

  script:
  """
  /local/two/Software/bitbucket/phylogeny/alignment_pruner.pl --file $aln --chi2_prune f${params.prune_percentage} > ${aln.baseName}.pruned.aln
  """
}

process buildTree  {
  input:
  file aln from pruned_alignment

  output:
  file "${aln.baseName}.treefile" into treefile
  file "RAxML_bipartitions.${aln.baseName}.PROTCATLG" optional true into raxml_tree

  cpus "${params.cpus}"
  publishDir "${params.outdir}", mode: 'copy'

  script:
  if (params.phylo_method == "raxml")
    """
    raxmlHPC-PTHREADS-SSE3 -s $aln -T ${task.cpus} -f a -x 12345 -p 12345 -N 100 -m PROTCATLG -n ${aln.baseName}.PROTCATLG
    cp "RAxML_bipartitions.${aln.baseName}.PROTCATLG" ${aln.baseName}.treefile
    """
  else if (params.phylo_method == "fasttree")
    """
    FastTree ${aln} > ${aln.baseName}.treefile
    """
  else if (params.phylo_method == "iqtree")
    """
    iqtree-omp -s ${aln} -m LG+C60 -bb 1000 -nt AUTO -pre ${aln.baseName}
    """
}
