
params.genomes = ""
params.references = ""
params.faas = ""
params.phylo_method = "fasttree"
params.outdir = "outdir"
params.numprots = 5
params.cpus = 10
params.prune_percentage = "0.0"

params.get_hits = true

optional_channel_boolean(params.get_hits).into{get_hits_database; get_hits_references}
additional_faas = optional_channel_from_path(params.faas)
optional_channel_from_path(params.references).into{references; references_add}
genomes = optional_channel_from_path(params.genomes)

def startup_message() {
    log.info "=========================================================="
    log.info "                       RP15 on Genomes"
    log.info "Author                            : Max Emil Schön"
    log.info "email                             : max-emil.schon@icm.uu.se"
    log.info "=========================================================="
    log.info "                       Parameters"
    log.info ""
    log.info "Genomes fasta file location       : $params.genomes"
    log.info "References COG location           : $params.references"
    log.info "Additional proteomes to add       : $params.faas"
    log.info "Output directory                  : $params.outdir"
    log.info "Phylogenetic method to use        : $params.phylo_method"
    log.info "Minimum number of proteins        : $params.numprots"
    log.info "number of CPUs                    : $params.cpus"
    log.info "heterogeneous sites to be removed : $params.prune_percentage"
    log.info "=========================================================="
    log.info "                      Usage examples:"
    log.info "only perform annnotation of the genomes:"
    log.info "nextflow run rp15_on_genomes.nf --genomes genomes/* --outdir outdir --get_hits false"
    log.info "=========================================================="
}

startup_message()


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

process cleanAdditionalFaa{
  input:
  file faa from additional_faas

  output:
  file "${faa.simpleName}.clean.faa" into additional_faas_clean

  tag {"${faa.simpleName}"}

  script:
  """
  #!/usr/bin/env anapy3

  from Bio import SeqIO
  with open("${faa.simpleName}.clean.faa", 'w') as outhandle:
      count = 0
      for rec in SeqIO.parse("$faa", 'fasta'):
          rec.id = "${faa.simpleName}_{}".format(count)
          count = count + 1
          SeqIO.write(rec, outhandle, 'fasta')
  """
}

faa_combined = faas.concat(additional_faas_clean).collectFile(name: 'all.faa', newLine: true)

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
  anapy3 $workflow.projectDir/concatenate.py *.aln -t $params.numprots
  """
}

process pruneAlignment {
  input:
  file aln from concatenated_alignment

  output:
  file "${aln.baseName}.pruned${params.prune_percentage}.aln" into pruned_alignment

  publishDir "${params.outdir}", mode: 'copy'

  script:
  """
  /local/two/Software/bitbucket/phylogeny/alignment_pruner.pl --file $aln --chi2_prune f${params.prune_percentage} > ${aln.baseName}.pruned${params.prune_percentage}.aln
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
    iqtree -s ${aln} -m LG+C60 -bb 1000 -nt AUTO -ntmax ${task.cpus} -pre ${aln.baseName}
    """
}

process makeNexus {
  input:
  file treefile

  output:
  file "${treefile.baseName}.nex" into nexus_file

  publishDir "${params.outdir}", mode: 'copy'

  script:
  """
  anapy3 $workflow.projectDir/create_nexus.py -i $treefile -o ${treefile.baseName}.nex
  """
}



def optional_channel_from_path(argument) {
  if(argument != ""){
    return Channel.fromPath(argument)
  }else{
    return Channel.empty()
  }
}

def optional_channel_boolean(argument) {
  if (argument){
    return Channel.from(Boolean.toString(argument))
  }else{
    return Channel.empty()
  }
}
