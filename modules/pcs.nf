include {strmem} from '../modules/fct_groovy.nf'

// We do PCA on qc2 data because relatively few SNPs and individuals will be removed later and
// this is an expensive operation so we start early. Similarly for computing relatedness
process compute_pcs {
  cpus params.max_plink_cores
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3
   input:
      path(plinks)
      val(pcsn)
      val(outputdir)
   output:
      tuple path("${prune}.eigenval"), path("${prune}.eigenvec"), emit:eigen
      tuple path ("${prune}.bed"), path("${prune}.bim"), path("${prune}.fam"), emit : plink_prune
   publishDir "${params.output_dir}/qc/pca", overwrite:true, mode:'copy',pattern: "${prune}*"
   script:
      base = plinks[0].baseName
      prune= "${base}-prune".replace(".","_")
      pcs = [pcsn,20].max()
     """
     plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
     plink --bfile ${base} --extract check.prune.in --make-bed --out $prune
     /bin/rm check*
     plink --bfile ${prune} --pca $pcs --out $prune
     """
}

process draw_pcs{
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3
   input:
      tuple path(eigvals), path(eigvecs)
      path(cc)
      val(col)
      val(diffpheno)
   output:
      tuple  path ("eigenvalue.pdf"), path(outputval)
   publishDir "${params.output_dir}/qc/pca/", overwrite:true, mode:'copy',pattern: "*.pdf"
   script:
      base=eigvals.baseName
      cc_fname = col
      // also relies on "col" defined above
      outputval="${base}-pca".replace(".","_")+".pdf"
      template "drawPCA.py"

}
