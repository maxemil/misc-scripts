params.alignments = ""
params.output_trees = ""
params.phylo_method = "iqtree"

clusters_alignments_trimmed = Channel.fromPath(params.alignments)

process cluster2Tree {
  input:
  file aln from clusters_alignments_trimmed

  output:
  file "${aln.simpleName}.*" into treefiles

  publishDir "${params.output_trees}", mode: 'copy'
  tag {"${aln.simpleName}"}
  cpus 3

  script:
  if (params.phylo_method == "iqtree")
    """
    iqtree -s $aln \
  	 -pre ${aln.simpleName} \
  	 -bb 1000 -bnni \
  	 -nt AUTO -ntmax ${task.cpus} \
  	 -m TESTNEW \
  	 -mset LG \
  	 -madd LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 \
  	 -keep-ident \
  	 -seed 12345 \
  	 -wbtl
    """
  else if (params.phylo_method == "phylobayes")
    """
    check_convergence(){
      while IFS= read -r ess
      do
        [ \$(echo \$ess">=100" | bc -l) -eq 1 ] || return 1
      done <<< "\$2"
      [ \$(echo \$1"<=0.3" | bc -l) -eq 1 ] || return 1
      return 0
    }

    mpirun -np 4 pb_mpi -d $aln -cat -lg "${aln.simpleName}".chain1 &
    pid1=\$!
    mpirun -np 4 pb_mpi -d $aln -cat -lg "${aln.simpleName}".chain2 &
    pid2=\$!

    sleep 60s

    while ps | grep -q \$pid1 && ps | grep -q \$pid1
    do
      burnin=\$(echo \$(tail -n 1 ${aln.simpleName}.chain2.trace | cut -f1)"/10" | bc)
      maxdiff=\$(bpcomp -x \$burnin ${aln.simpleName}.chain2.treelist ${aln.simpleName}.chain1.treelist | grep maxdiff | cut -d ':'  -f 2)
      tracecomp_ess=\$(tracecomp -x \$burnin ${aln.simpleName}.chain1.trace ${aln.simpleName}.chain2.trace | sed "s/[\\t ]\\+/\\t/g"  | cut -f2 | grep -o '[0-9]\\+')
      if check_convergence \$maxdiff \$tracecomp_ess
      then
        echo 0 > ${aln.simpleName}.chain1.run
        echo 0 > ${aln.simpleName}.chain2.run
        echo \$(date) \$maxdiff "chains converged, stopping" >> ${aln.simpleName}.bpcomp.log
      else
        echo \$(date) \$maxdiff >> ${aln.simpleName}.bpcomp.log
      fi
      sleep 60m
    done
    """
}
