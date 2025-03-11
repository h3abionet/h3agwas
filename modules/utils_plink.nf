include {strmem} from '../modules/fct_groovy.nf'
process clean_plink{
 cpus params.max_cpus
 input:
  tuple path(bed), path(bim), path(fam)
  path(snp_exclude)
  path(snp_include)
  path(indkeep)
  val(maf)
  val(plink_thin_count)
  val(just_acgt)
 output:
  tuple path("${outfile}.bed"),  path("${outfile}.bim"), path("${outfile}.fam")
 script:
  bfile=bed.baseName
  outfile=params.output+"_clean"
  snp_exclude=snp_exclude.toString()
  snp_include=snp_include.toString()
  range=(snp_exclude=="01" || snp_exclude=="02" || snp_exclude=="03") ? "" : " --exclude range $snp_exclude "
  range2=(snp_include=="01" || snp_include=="02" ||snp_include=="03") ? "" : " --extract range $snp_include "
  mafcmd=(maf==0) ? "" : " --maf ${maf}"
  cleanrsoption=(just_acgt==1) ?  " --snps-only just-acgt " : ""
  """
  plink --bfile $bfile --threads ${params.max_cpus}  $range --make-bed --out $outfile $mafcmd --keep $indkeep --keep-allele-order $range2 $cleanrsoption --allow-extra-chr
  """
}

process getfrequency{
 input :
     tuple path(bed), path(bim), path(fam)
 publishDir "${params.output_dir}/checkmichigan/frq", overwrite:true, mode:'copy'
 output :
     tuple path(bim),path("${bfile}.frq")
 script :
    bfile=bed.baseName
    """
    plink --freq -bfile $bfile  -out $bfile
    """
 }

 process list_chro {
  input :
      tuple path(bed), path(bim), path(fam)
  publishDir "${outputdir}/", overwrite:true, mode:'copy'
  output :
     stdout
  script :
      """
      awk '{print \$1}' $bim|uniq|sort|uniq
      """
 }

process indep_pairwise {
  cpus params.max_cpus
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3
   input:
      path(plinks)
      val(indeppairwise)
      val(outputdir)
   output:
      tuple path ("${prune}.bed"), path("${prune}.bim"), path("${prune}.fam")
   publishDir "${outputdir}/", overwrite:true, mode:'copy',pattern: "${prune}*"
   script:
      base = plinks[0].baseName
      prune= "${base}-prune".replace(".","_")
     """
     plink --bfile ${base} --indep-pairwise ${indeppairwise} --out check 
     plink --bfile ${base} --extract check.prune.in --make-bed --out $prune 
     """
} 
