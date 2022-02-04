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
 * Description : pipeline to do a Conditional analysis using gemnma
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths

def helps = [ 'help' : 'help' ]
allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "AMI", "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "other_mem_req", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin"]

params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.output = "out_cond"
params.thin       = ""
params.covariates = ""
params.chro_cond      = -1
params.pos_cond      = -1
params.around = 250000
params.pos = ""
params.pos_ref = -1
params.file_rs_buildrelat = ""
params.rs_list=""
params.sample_snps_rel=0
params.gemma_num_cores=4
params.gemma_mem_req="4GB"
params.gemma_bin = "gemma"
params.gemma_mat_rel=""
params.gemma_relopt = 1
params.max_plink_cores=5
params.plink_mem_req = '6GB' // how much plink needs for this
params.sample_snps_rel_paramplkl="100 20 0.1 --maf 0.01"




max_plink_cores = params.max_plink_cores
plink_mem_req = params.plink_mem_req

filescript=file(workflow.scriptFile)
projectdir="${filescript.getParent()}"
dummy_dir="${projectdir}/../qc/input"


def fileColExists = { fname, pname, cname ->
  if (fname.contains("s3://")){
    println "The file <$fname> is in S3 so we cannot do a pre-check";
    return;
  }
  if (fname.contains("az://")){
    println "The file <$fname> is in Azure so we cannot do a pre-check";
    return;
  }
  f = new File(fname)
  if (! f.exists() && ! fname.contains("s3://") && ! fname.contains("az://")) {
     error("\n\nThe file <${fname}> given for <${pname}> does not exist")
    } else {
      def line
      f.withReader { line = it.readLine() }
      // now get the column headers
      fields = line.split()
      // now separate the column
      cols = cname.split(",")
      cols.each { col ->
        det = col.split("/")
        if ((det[0].length()>0) && (! fields.contains(det[0])))
          error("\n\nThe file <${fname}> given for <$pname> does not have a column <${det}>\n")
      }
    }
}

fileColExists(params.data,"${params.data} - covariates", params.covariates)
fileColExists(params.data,"${params.data} - phenotypes", params.pheno)

//---- Modification of variables for pipeline -------------------------------//


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



bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")


Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }
extract_reg_ch=Channel.create()
rel_ch_gemma=Channel.create()
ch_select_rs_format=Channel.create()
raw_src_ch.separate( extract_reg_ch,rel_ch_gemma , ch_select_rs_format) { a -> [a,a,a] }



process extract_region{
    input:
      set file(bed), file(bim), file(fam) from  extract_reg_ch
    output:
      /*JT Append initialisation boltlmm_assoc_ch */
      set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into  (gemma_assoc_ch, pheno_assoc_ch, gemmacond_assoc_ch, ld_assoc_ch, gemmacondall_assoc_ch)
    script:
       base = bed.baseName
       out  = base+"_t"
       chroC=params.chro_cond
       posCmin=params.pos_cond.toString().split(',').collect{it as int }.min()
       posCmax=params.pos_cond.toString().split(',').collect{it as int }.max()
       posCmin=(params.pos_ref<0) ? "$posCmin" : "${[posCmin, params.pos_ref].min()}"
       posCmax=(params.pos_ref<0) ? "$posCmax" : "${[posCmax, params.pos_ref].max()}"
       filerange="tmp.ped"
       """
       begin=`expr $posCmin-${params.around}`
       end=`expr $posCmax + ${params.around}`
       if [ "\$begin" -lt 1 ]
       then
       begin=1
       fi
       echo -e "$chroC\t\$begin\t\$end\t$chroC:$posCmin" > $filerange
       plink --keep-allele-order --bfile $base --make-bed --out $out --extract range $filerange
       """
}

pheno_file_ch = Channel.fromPath(params.data, checkIfExists:true)
process add_condpheno{
    label 'R'
    input:
      set file(bed), file(bim), file(fam) from pheno_assoc_ch
      file(filepheno) from pheno_file_ch
    output :
      file(pheno) into pheno_filecond_form_ch, pheno_file_form_ch,pheno_filecondall_form_ch
      stdout into (rsname_ch, rsname_ch2)
    script :
        pheno=filepheno+'_condadd.pheno'
        bfile=bed.baseName
        """
        cond_changepheno.r --bfile $bfile --chro_cond ${params.chro_cond} --pos_cond $params.pos_cond --data $filepheno --out $pheno &> Rvalue
        awk '{print \$2}' tmp.map
        """
}
check=Channel.create()
rsname_ch_f=Channel.create()
rsname_ch.flatMap { list_str -> list_str.split() }.tap ( check) .set {rsname_ch_f}





if(params.gemma_mat_rel==""){
    if(params.file_rs_buildrelat=="" && params.sample_snps_rel==1){
      process select_rs_format{
       cpus max_plink_cores
       memory plink_mem_req
       time   params.big_time
       input :
         set file(bed),file(bim), file(fam) from ch_select_rs_format
      output:
       file("${prune}.prune.in") into   filers_matrel_mat_gem
      script:
        base = bed.baseName
        prune= "${base}-prune"
        """
        plink --bfile ${base} --indep-pairwise ${params.sample_snps_rel_paramplkl} --out $prune   --threads ${params.max_plink_cores}
        """
    }
    balise_filers_rel=1
   }else{
     if(params.file_rs_buildrelat==""){
        filers_matrel_mat_gem=file("${dummy_dir}/0") 
        balise_filers_rel=0
     }else{
        filers_matrel_mat_gem=Channel.fromPath(params.file_rs_buildrelat)
        balise_filers_rel=1
     }
   }
    process getGemmaRel {
      label 'gemma'
      cpus params.gemma_num_cores
      memory params.gemma_mem_req
      time params.big_time
      input:
        file plinks from rel_ch_gemma
        file file_rs from filers_matrel_mat_gem
      publishDir "${params.output_dir}/gemma/rel", overwrite:true, mode:'copy'
      output:
         file("output/${base}.*XX.txt") into rel_mat_ch, rel_matcondall_ch, rel_matcond_ch
      script:
        base = plinks[0].baseName
        famfile=base+".fam"
        rs_list = balise_filers_rel==1 ? " -snps $file_rs " : ""
        """
        export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
        cat $famfile |awk '{print \$1"\t"\$2"\t"0.2}' > pheno
        ${params.gemma_bin} -bfile $base  -gk ${params.gemma_relopt} -o $base -p pheno -n 3 $rs_list
       """
    }  


  }else{
   rel_matcond_ch=Channel.fromPath(params.gemma_mat_rel)
   rel_matcondall_ch=Channel.fromPath(params.gemma_mat_rel)
   rel_mat_ch=Channel.fromPath(params.gemma_mat_rel)
  }
   if(params.rs_list==""){
        rsfile=file("${dummy_dir}/0")
        rsfilecond=file("${dummy_dir}/0")
        rsfilecondall=file("${dummy_dir}/0")
     }else{
        rsfile=file(params.rs_list)
        rsfilecond=file(params.rs_list)
        rsfilecondall=file(params.rs_list)
   }

 process doGemmaCond {
    maxForks params.max_forks
    label 'gemma'
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      file(pheno) from pheno_filecond_form_ch
      file(rel) from rel_matcond_ch
      file(plinks) from  gemmacond_assoc_ch
      file(rsfilelist) from rsfilecond
    each rsname from rsname_ch_f
    publishDir "${params.output_dir}/$out", overwrite:true, mode:'copy'
    output:
      file("${dir_gemma}/${out}.log.txt")
      file("${dir_gemma}/${out}.assoc.txt") into gemma_outcond_ch
    script:
       this_pheno         = "${params.pheno}"
       our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
       our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
       our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
       data_nomissing     = "pheno-"+our_pheno+".pheno"
       list_ind_nomissing = "lind-"+our_pheno+".lind"
       rel_matrix         = "newrel-"+our_pheno+".rel"
       base               =  plinks[0].baseName
       inp_fam            =  base+".fam"
       newbase            =  base+"-"+our_pheno
       newfam             =  newbase+".fam"
       gemma_covariate    = "${newbase}.gemma_cov"
       phef               = "${newbase}_n.phe"
       covariate_option = " --cov_list $rsname"
       covariate_option = (params.covariates=="") ? " $covariate_option " : "$covariate_option,${params.covariates}" 
       covar_opt_gemma    =  " -c $gemma_covariate " 
       rs_plk_gem         =  (params.rs_list) ?  " --extract  $rsfilelist" : ""
       out                = rsname.replace(':','_')
       dir_gemma          =  "gemma"
       """
       list_ind_nomissing.py --data $pheno --inp_fam $inp_fam $covariate_option --pheno $our_pheno3 --dataout $data_nomissing \
                             --lindout $list_ind_nomissing
       gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
       plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase  ${rs_plk_gem}
       all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
                          --pheno $our_pheno2 --phe_out ${phef} --form_out 1
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf 0.0000001
       mv output ${dir_gemma}
       rm $rel_matrix
       rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
       """
  }

 process doGemmaCondAllSnp {
    maxForks params.max_forks
    label 'gemma'
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      file(pheno) from pheno_filecondall_form_ch
      file(rel) from rel_matcondall_ch
      file(plinks) from  gemmacondall_assoc_ch
      file(rsfilelist) from rsfilecondall
      val(listrs) from rsname_ch2
    publishDir "${params.output_dir}/$out", overwrite:true, mode:'copy'
    output:
      file("${dir_gemma}/${out}.log.txt")
      file("${dir_gemma}/${out}.assoc.txt") into gemma_outcondall_ch
    script:
       this_pheno         = "${params.pheno}"
       our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
       our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
       our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
       data_nomissing     = "pheno-"+our_pheno+".pheno"
       list_ind_nomissing = "lind-"+our_pheno+".lind"
       rel_matrix         = "newrel-"+our_pheno+".rel"
       base               =  plinks[0].baseName
       inp_fam            =  base+".fam"
       newbase            =  base+"-"+our_pheno
       newfam             =  newbase+".fam"
       gemma_covariate    = "${newbase}.gemma_cov"
       phef               = "${newbase}_n.phe"
       listrs=listrs.replace(' ',',').replace('\n',',').replaceAll(/,$/,'')
       covariate_option = " --cov_list $listrs"
       covariate_option = (params.covariates=="") ? " $covariate_option " : "$covariate_option,${params.covariates}"
       covar_opt_gemma    =  " -c $gemma_covariate "
       rs_plk_gem         =  (params.rs_list) ?  " --extract  $rsfilelist" : ""
       out                = "allrs"
       dir_gemma          =  "gemma"
       println listrs
             """
       list_ind_nomissing.py --data $pheno --inp_fam $inp_fam $covariate_option --pheno $our_pheno3 --dataout $data_nomissing \
                             --lindout $list_ind_nomissing
       gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
       plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase  ${rs_plk_gem}
       all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
                          --pheno $our_pheno2 --phe_out ${phef} --form_out 1
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf 0.0000001
       mv output ${dir_gemma}
       rm $rel_matrix
       rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
       """
  }

  

 process doGemma {
    maxForks params.max_forks
    label 'gemma'
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      file(pheno) from pheno_file_form_ch
      file(rel) from rel_mat_ch
      file(plinks) from  gemma_assoc_ch
      file(rsfilelist) from rsfile
    publishDir "${params.output_dir}/initial/", overwrite:true, mode:'copy'
    output:
      file("${dir_gemma}/${out}.log.txt")
      file("${dir_gemma}/${out}.assoc.txt") into gemma_out_ch
    script:
       this_pheno         = "${params.pheno}"
       our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
       our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
       our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
       data_nomissing     = "pheno-"+our_pheno+".pheno"
       list_ind_nomissing = "lind-"+our_pheno+".lind"
       rel_matrix         = "newrel-"+our_pheno+".rel"
       base               =  plinks[0].baseName
       inp_fam            =  base+".fam"
       newbase            =  base+"-"+our_pheno
       newfam             =  newbase+".fam"
       gemma_covariate    = "${newbase}.gemma_cov"
       phef               = "${newbase}_n.phe"
       covariate_option = (params.covariates=="") ? " " : "  --cov_list  ${params.covariates}"
       covar_opt_gemma    =  (params.covariates=="") ? " " :" -c $gemma_covariate "
       rs_plk_gem         =  (params.rs_list) ?  " --extract  $rsfilelist" : ""
       out                = 'intial'
       dir_gemma          =  "gemma"
       """
       list_ind_nomissing.py --data $pheno --inp_fam $inp_fam $covariate_option --pheno $our_pheno3 --dataout $data_nomissing \
                             --lindout $list_ind_nomissing
       gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
       plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase  ${rs_plk_gem}
       all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
                          --pheno $our_pheno2 --phe_out ${phef} --form_out 1
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf 0.0000001
       mv output ${dir_gemma}
       rm $rel_matrix
       rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
       """
  }

process computed_ld{
    input :
      file(plinks) from ld_assoc_ch
    output :
      file(ld) into ld_res_notsq
      file(ld2) into ld_res_sq
      set file("${newbase}.bim"), file("${newbase}.bed"),file("${newbase}.fam")  into plk_pos_ch
    publishDir "${params.output_dir}/ld/", overwrite:true, mode:'copy'
    script :
       listpos= (params.pos_ref<1) ?"${params.pos_cond} "  : "${params.pos_cond},${params.pos_ref}"
       chr=params.chro_cond
       base               =  plinks[0].baseName
       out = params.output
       ld=params.output+'_notsq.ld'
       ld2=params.output+'_sq.ld'
       newbase=params.output+'_sub'
       """
       echo $listpos|awk -v chro=$chr -F"," '{for(Cmt=1;Cmt<=NF;Cmt++)print chro"\t"\$Cmt"\t"\$Cmt"\t"chro":"\$Cmt}' >ldlist 
       plink -bfile $base --make-bed --extract range ldlist -out $newbase --keep-allele-order
       plink -bfile $newbase  --r2           --ld-window 99999  --out ${params.output}_notsq    --ld-window-kb 10000 --ld-window-r2 0
       plink -bfile $newbase  --r2 square --out ${params.output}_sq     
       """
}

process plot_ld{
   label 'R'
  input : 
     file(ld2) from ld_res_sq
     set file(bim), file(bed), file(fam) from plk_pos_ch
  output :
     file(fileout) into plot_ld
  publishDir "${params.output_dir}/res/fig/", overwrite:true, mode:'copy'
  script :
      fileout=params.output+'_ld.pdf'
      posref= (params.pos_ref<1) ?""  : " --pos_ref  ${params.pos_ref}"
      """
      cond_plotld.r --ld $ld2 --bim $bim --out $fileout  $posref
      """
}

combine_gemm_ch=gemma_outcond_ch.collect()
process plot_res{
   label 'R'
   input :
      file(ld) from ld_res_notsq
      file(gemmai) from gemma_out_ch
      file(gemmamerge) from gemma_outcondall_ch
      file(allgemmars) from combine_gemm_ch
   publishDir "${params.output_dir}/res/fig/", overwrite:true, mode:'copy'
   output :
      file("${out}*")
   script :
     allrsgemma=allgemmars.join(',')
     out=params.output
     """
      cond_plotresgemma.r --ld $ld --list_file_assoc $allrsgemma --file_I $gemmai --file_merge $gemmamerge --pos_ref $params.pos_ref  --out $out
     """ 
}


