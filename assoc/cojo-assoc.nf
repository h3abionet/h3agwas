#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2019
 *
 *
 * Description : pipeline to do a Conditional and joint multiple-SNP analysis of GWAS
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths



def helps = [ 'help' : 'help' ]

allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "AMI", "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "other_mem_req", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin"]
// define param for
//annotation_model=["gemma","boltlmm", "plink", "head", "linear", "logistic", "fisher", "fastlmm", ""]
//allowed_params+=annotation_model
annotation_param=[ "file_gwas","gcta_bin","cojo_p","cojo_p","gcta_mem_req","cojo_slct_other", "cojo_top_snps","cojo_slct", "cojo_actual_geno"]
allowed_params+=annotation_param
allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2", "head_n"]
allowed_params+=allowed_params_head


params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}

def params_help = new LinkedHashMap(helps)


params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.output_testing = "cleaned"
params.covariates = ""
params.gcta_bin="gcta64"
outfname = params.output_testing

params.head_pval = "P_BOLT_LMM"
params.head_freq = ""
params.head_n = ""
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.head_beta="BETA"
params.head_se="SE"
params.head_A1="ALLELE0"
params.head_A2="ALLELE1"
params.cojo_p=1e-7
params.cojo_wind=10000
params.cut_maf=0.01
params.gcta_mem_req="15GB"
params.plink_mem_req="6GB"
params.gcta_cpus_req = 1
params.cojo_slct=1
params.cojo_slct_other=""
params.cojo_actual_geno=0
params.cojo_top_snps=0
params.big_time='100h'



gcta_mem_req=params.gcta_mem_req
gcta_cpus_req = params.gcta_cpus_req+2
plink_mem_req = params.plink_mem_req
max_plink_cores = params.max_plink_cores
params.help = false

def getConfig = {
  all_files = workflow.configFiles.unique()
  text = ""
  all_files.each { fname ->
      base = fname.baseName
      curr = "\n\n*-subsection{*-protect*-url{$base}}@.@@.@*-footnotesize@.@*-begin{verbatim}"
      file(fname).eachLine { String line ->
        if (line.contains("secretKey")) { line = "secretKey='*******'" }
        if (line.contains("accessKey")) { line = "accessKey='*******'" }
        curr = curr + "@.@"+line
      }
      curr = curr +"@.@*-end{verbatim}\n"
      text = text+curr
  }
  return text
}

// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}

println "\nTesting data : ${params.data}\n"
println "Testing for gwas file : ${params.file_gwas}\n"



bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()

raw_src_ch= Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }

plink_format=Channel.create()
plink_cojo=Channel.create()
plink_other=Channel.create()
raw_src_ch.separate (plink_format, plink_cojo, plink_other) { a -> [ a, a, a] }


gwas_chrolist = Channel.fromPath(params.file_gwas)
//COunt list chro
process getListeChro{
        input :
          file(gwas_res) from gwas_chrolist
        output :
          stdout into (chrolist,chrolist2)
        script:
         """
         PosChro=`head -1 ${gwas_res}|awk \'{for(Cmt=1;Cmt<=NF;Cmt++)print Cmt\"\\t\"\$Cmt}\'|grep -e \"\t\"${params.head_chr}\$|awk \'{print \$1}\'`
         sed \'1d\' ${gwas_res}|awk -v poschro=\$PosChro \'{print \$poschro}\'|uniq|sort|uniq
        """
}
check2 = Channel.create()
chrolist=chrolist.flatMap { list_str -> list_str.split() }.tap ( check2)

if(params.data){
data_listind=Channel.fromPath(params.data)
process getListInd{
   input :
     file data from data_listind
   output:
     file(keepout) into (filekeepformat,filekeepcojo)
   script :
     varpheno=params.pheno!=""? " --pheno ${params.pheno}" : ""
   keepout="list_ind.keep"
   """
   gcta_extindplk.py --data ${data}  $varpheno --out $keepout
   """
}
}else{
famlind = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()
process getListInd2{
   input :
     file(fam) from famlind
   output :
     file(keepout) into (filekeepformat,filekeepcojo)
   script :
   keepout="list_ind.keep"
   """
   awk \'{print \$1\" \"\$2}\' $fam > $keepout
   """

}
//filekeepcojo=filekeepformat=file('NO_FILEINDKEEP')
}


gwas_format = Channel.fromPath(params.file_gwas)
chrolist_ch=Channel.create()
check = Channel.create()
chrolist.flatMap { list_str -> list_str.split() }.tap ( check) .set {chrolist_ch}


process doFormatData{
   cpus params.max_plink_cores
   memory params.gcta_mem_req
   input :
     file(keepind) from filekeepformat
     set file(bed), file(bim), file(fam) from plink_format
     file(gwas) from gwas_format
   each chro from chrolist_ch 
   output :
        set file(keepind), file(bed), file(bim), file(fam),chro, file(out) into gwas_chro_cojo,gwas_chro_topsnp
   script :
      out=chro+".format"
      baseplk=bed.baseName
      headfreq=params.head_freq!="" ? " --freq_header ${params.head_freq}" : ""
      headn=params.head_n!="" ? " --freq_header ${params.head_n}" : ""
      headkeep=params.data!="" ? " --keep $keepind " : ""
      """
      gcta_format.py --inp_asso $gwas  --rs_header ${params.head_rs} --pval_header ${params.head_pval} $headfreq --a1_header ${params.head_A1} --a2_header ${params.head_A2} --se_header ${params.head_se} --beta_header ${params.head_beta} --chro_header ${params.head_chr} --chr $chro --bfile $baseplk --out $out --threads ${params.max_plink_cores}
      """
}
if(params.cojo_slct){
process SLCTAnalyse{
    time params.big_time
    //memory { gcta_mem_req * task.attempt }
    memory gcta_mem_req
    cpus gcta_cpus_req
    //errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    //memory gcta_mem_req
    input :
        set file(keepind), file(bed), file(bim), file(fam),chro, file(gwas_chro) from gwas_chro_cojo
    output :
        file("${out}.jma.cojo")  optional true into res_chro_slct
    script :
      baseplk=bed.baseName
      out="res_cojo"+chro
      headkeep=params.data!="" ? " --keep $keepind " : ""
      cojoactu=params.cojo_actual_geno==1 ? " --cojo-actual-geno " : ""
      """
      export OMP_NUM_THREADS=${params.gcta_cpus_req}
      ${params.gcta_bin} --bfile $baseplk --chr $chro --out $out --cojo-file ${gwas_chro} --cojo-p ${params.cojo_p} --maf ${params.cut_maf} ${headkeep} --cojo-slct ${params.cojo_slct_other} --cojo-wind ${params.cojo_wind} $cojoactu --thread-num ${params.gcta_cpus_req} &> $out".out"
      gcta_manage_error.py --file_err $out".out"
      if [ ! -f ${out}.jma.cojo ]
      then 
      echo -e "Chr\\tSNP\\tbp\\trefA\\tfreq\\tb\\tse\\tp\\tn\\tfreq_geno\\tbJ\\tbJ_se\\tpJ\\tLD_r" > ${out}.jma.cojo
      fi
      """ 
}
res_chro_slct_merg=res_chro_slct.collect()
process SLCTMerge{
   input :
      file(list_file) from res_chro_slct_merg
   publishDir "${params.output_dir}/slct", overwrite:true, mode:'copy'
   output :
      file("$out") into res_slct
   script : 
     fnames = list_file.join(" ")
     file1  = list_file[0]
     out=params.output+"_slct.jma.out"
     """
     head -1 $file1 > $out
     tail -q -n +2  $fnames >> $out
     """
}
//to do report
}

if(params.cojo_top_snps_chro>0){
process TopAnalyse{
    time  params.big_time
    //memory { gcta_mem_req * task.attempt }
    memory gcta_mem_req
    cpus gcta_cpus_req
    //errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    input :
        set file(keepind), file(bed), file(bim), file(fam),chro, file(gwas_chro) from gwas_chro_topsnp
    output :
        file("${out}.jma.cojo") into res_chro_top
    script :
      baseplk=bed.baseName
      out="res_cojo"+chro+"."+params.cojo_top_snps_chro
      headkeep=params.data!="" ? " --keep $keepind " : ""
      cojoactu=params.cojo_actual_geno==1 ? " --cojo-actual-geno " : ""
      """
      export OMP_NUM_THREADS=${params.gcta_cpus_req}
      ${params.gcta_bin} --bfile $baseplk --chr $chro --out $out --cojo-file ${gwas_chro} --maf ${params.cut_maf} ${headkeep} --cojo-top-SNPs $params.cojo_top_snps_chro ${params.cojo_slct_other} --cojo-wind ${params.cojo_wind} $cojoactu --thread-num ${params.gcta_cpus_req}
      """
}
res_chro_top_merg=res_chro_top.collect()
process TopMerge{
   input :
      file(list_file) from res_chro_top_merg
   publishDir "${params.output_dir}/top", overwrite:true, mode:'copy'
   output :
      file("$out") into res_top
   script :
     fnames = list_file.join(" ")
     file1  = list_file[0]
     out=params.output+"_top.jma.out"
     """
     head -1 $file1 > $out
     tail -q -n +2  $fnames >> $out
     """
}
//to do report
}









