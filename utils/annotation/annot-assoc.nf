#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Shaun Aron
 *   	Rob Clucas
 *      Eugene de Beste
 *      Lerato Magosi
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description : pipeline annotation 
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths




def helps = [ 'help' : 'help' ]

allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "AMI", "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "other_mem_req", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin"]
// define param for 
annotation_model=["gemma","boltlmm", "plink", "head", "linear", "logistic", "fisher", "fastlmm", ""]
allowed_params+=annotation_model
annotation_param=["list_rs", "file_annotation", "file_gwas", "around_rs"]
allowed_params+=annotation_param
allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2"]
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
outfname = params.output_testing


/* Defines the path where any scripts to be executed can be found.
 */

/**/
params.head_pval = "P_BOLT_LMM"
params.head_freq = ""
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.head_beta="BETA"
params.head_se="SE"
params.head_A1="ALLELE1"
params.head_A2="ALLELE0"

params.around_rs=100000
params.cut_maf = 0.01

params.loczm_bin  = ""
params.loczm_pop = "AFR"
params.loczm_build = "hg19"
params.loczm_source ="1000G_Nov2014"
params.loczm_gwascat = ""
params.data=""
params.pheno=""
params.list_rs=""

params.list_file_annot=""
params.info_file_annot=""

if(params.list_rs==""){
error("params : list_rs not initialise, must be rs contains in bed / gwas file")
}
if(params.data==""){
error("params : data not initialise, must be data file ")
}
if(params.pheno==""){
error("params : pheno not initialise, must be column of data file")
}



/*JT Append initialisation variable*/
params.bolt_impute2fidiid=""
/*gxe param : contains column of gxe*/
params.gxe=""

params.input_pat  = 'raw-GWA-data'
//format boltlmm
params.boltlmm = 0

params.sexinfo_available = "false"


params.plink_mem_req = '750MB' // how much plink needs for this
params.other_process_memory = '8GB' // how much other processed need


plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_memory
max_plink_cores = params.max_plink_cores 
params.help = false
//ListeDB=(refGene cg46 cg69 dbnsfp31a_interpro dbscsnv11 esp6500siv2_all exac03 gme intervar_20170202 knownGene mcap avsnp150 ljb26_all dbnsfp33a 1000g2015aug kaviar_20150923 avsnp150 clinvar_20170905 LoFtool)
params.list_db_anovar="1000g2015aug,gme"
params.ref_annovar="hg19"
params.bin_annovar="annotate_variation.pl"


//data_ch = Channel.fromPath(params.data, checkIfExists:true)

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
if(params.boltlmm==1){
head_pval = "P_BOLT_LMM"
head_freq = "A1FREQ"
head_bp = "BP"
head_chr = "CHR"
head_rs = "SNP"
head_beta="BETA"
head_se="SE"
head_A1="ALLELE1"
head_A2="ALLELE0"

}else{
head_pval=params.head_pval 
head_freq=params.head_freq
head_bp=params.head_bp
head_chr=params.head_chr
head_rs=params.head_rs
head_beta=params.head_beta
head_se=params.head_se
head_A1=params.head_A1
head_A2=params.head_A2
}

if(params.loczm_bin==""){
error("\n\n------params loczm_bin must be initialise \n\n\n---\n")
}




println "\nTesting rs : ${params.list_rs}\n"
println "\nTesting data : ${params.data}\n"
println "Testing for phenotypes  : ${params.pheno}\n"
println "Testing for gwas file : ${params.file_gwas}\n"
println "Using covariates        : ${params.covariates}\n"
println "Around rs : ${params.around_rs}\n\n"

gwas_ch = Channel.fromPath(params.file_gwas, checkIfExists:true)

list_rs=params.list_rs.split(",")
rs_label_ch=Channel.from(list_rs)


process ExtractInfoRs{
    memory other_mem_req
    input:
       file(gwas_file) from gwas_ch
    each rs from rs_label_ch
    output :
      set val(rs), file(out_locus_rs) into locuszoom_ch
      set val(rs), file(out_info_rs) into infors_rs, infors_rs2
      set val(rs), file(out_gwas_rs) into infors_gwas
    script:
      out="sub_"+rs.replace(':',"_")
      out_locus_rs=out+"_around.stat"
      out_gwas_rs=out+"_gwas.stat"
      out_info_rs=out+"_info.stat"
      infoaf=head_freq=="" ? "" :  " --freq_header  $head_freq --maf ${params.cut_maf} " 
      """
     an_extract_rs.py --inp_resgwas  $gwas_file --chro_header $head_chr --pos_header $head_bp --rs_header $head_rs --pval_header $head_pval --beta_header ${head_beta} --list_rs $rs --around_rs ${params.around_rs} --out_head $out
      """
}

if(params.list_file_annot!=""){
fileannot_ch = infors_rs.join(Channel.from(list_rs).combine(Channel.fromPath(params.list_file_annot,checkIfExists:true))).join(Channel.from(list_rs).combine(Channel.fromPath(params.info_file_annot,checkIfExists:true)))
}else{

gwas_ch_annovar = Channel.fromPath(params.file_gwas,checkIfExists:true)
process ExtractAllRs{
    memory other_mem_req
    input:
       file(gwas_file) from gwas_ch_annovar
    output :
      file("${out}_all.annov") into info_annovar, info_annovar_chro, info_annovar_2
    script :
       out=params.output
       """ 
      an_extract_rs.py --inp_resgwas  $gwas_file --chro_header $head_chr --pos_header $head_bp --rs_header $head_rs --pval_header $head_pval --beta_header ${head_beta} --list_rs $params.list_rs --around_rs ${params.around_rs} --out_head $out --a1_header $head_A1 --a2_header  $head_A2
       """
}
  process getListeChro{
        input :
          file(BimFile) from info_annovar_chro
        output :
          stdout into chrolist
        script:
         """
         cat $BimFile|awk '{print \$1}'|uniq|sort|uniq
        """
  }
  chrolist=chrolist.flatMap { list_str -> list_str.split() }


listdb=Channel.from(params.list_db_anovar.split(','))
process getannovar_db {
  input : 
      val(db) from listdb
  publishDir "${params.output_dir}/dbzip/", overwrite:true, mode:'copy'
  output:
    file("humandb/${params.ref_annovar}_*.txt") into db_annovar
  script :
   output=params.ref_annovar+"_" + db+'.txt'
   """
   ${params.bin_annovar} -downdb -buildver ${params.ref_annovar} -webfrom annovar $db humandb
   """
}

info_annovar=info_annovar.combine(db_annovar.flatMap())

process annovar_chro_db{
   input :
     set file(filepos_annov), file(annov_db) from info_annovar
   each chro from chrolist
   output :
      set val(chro),file("$fileout") into list_db_format_ch_chro
   script :
     db=annov_db.baseName
     out_chro_tmp=db
     out_annov_tmp=chro+"_tmp.annov"
     fileout=db+"_"+chro+".txt"
     """
     awk -v chro=$chro '{if(\$1==chro)print \$0}' $filepos_annov > $out_annov_tmp
     head -1 $annov_db > $out_chro_tmp
     awk -v chro=$chro '{if(\$1==chro)print \$0}' $annov_db >> $out_chro_tmp
     an_extract_info.py --list_pos $out_annov_tmp --annov_file $out_chro_tmp --out $fileout 
     rm $out_chro_tmp $out_annov_tmp
     """
}

list_db_chro=info_annovar_2.combine(list_db_format_ch_chro.groupTuple())
process mergechro_db{
 input :
    set file(infoannov),val(chro), file(listfile) from list_db_chro
  publishDir "${params.output_dir}/db_annov/", overwrite:true, mode:'copy'
  output :
    file(out) into merge_chro_ch
  script :
    out_annov_tmp=chro+"_tmp.bed"
    listfch=listfile.join(",")
    out=chro+"_merge.annov"
    """
    awk -v chro=$chro '{if(\$1==chro)print \$0}' $infoannov > $out_annov_tmp
    merge_file.py --list_pos $out_annov_tmp --list_file $listfch --out $out
    """
}


db_format_ch=merge_chro_ch.collect()
process merge_annovar{
   input :
    file(listfile) from db_format_ch
   output:
    set file(file_annot), file(info_file_annot) into fileannot_ch  
   script :
      file_annot="annot"
      info_file_annot="annot.info"
      """
      ls -l ${params.output_dir}/db_annov/*_merge.annov | awk -F"[_/]" '{\${NF-1}"\t"\$0}'> $file_annot  
      touch $info_file_annot
      """
}
//fileannot_ch = infors_rs.join(Channel.from(list_rs).combine(Channel.fromPath(params.list_file_annot))).join(Channel.from(list_rs).combine(Channel.fromPath(params.info_file_annot)))
}


if(params.loczm_gwascat!=""){
loczm_gwascat=" --gwas-cat ${params.loczm_gwascat}"
}else{
loczm_gwascat=""
}

loczm_dir=file(params.loczm_bin).getParent().getParent()

locuszoom_ch=locuszoom_ch.combine(Channel.fromPath(loczm_dir,type:'dir'))
process PlotLocusZoom{
    label 'py2R'
    memory plink_mem_req
    input : 
       set val(rs), file(filegwas),file(lz_dir) from locuszoom_ch
    publishDir "${params.output_dir}/$rsnameout", overwrite:true, mode:'copy'
    output :
       file("out_*/*.svg")
       set val(rs), file("out_*/*.pdf") into report_lz_ch
    script :
       rsnameout=rs.replace(':',"_")
       """
       rs2=`echo $rs | awk -F':' '{if(NF==1){print \$0;}else{print \$1\":\"\$2;}}'` 
       ${lz_dir}/bin/locuszoom --epacts  $filegwas --delim tab --refsnp  \$rs2 --flank ${params.around_rs} --pop ${params.loczm_pop} --build ${params.loczm_build} --source ${params.loczm_source} $loczm_gwascat --svg  -p out --no-date 
       """
}

process ExtractAnnotation{
      label 'latex'
      memory plink_mem_req
      input :
        set val(rs),file(file_rs),file(annot_file), file(annot_info) from fileannot_ch
      publishDir "${params.output_dir}/$rsnameout", overwrite:true, mode:'copy'
      output :
        file("${out}*")
        file("*${rsnameout}*")
        set val(rs), file("${out}.pdf") into report_info_rs
      script :    
         out="annot-"+rs.replace(':','_')
         rsnameout=rs.replace(':',"_")
         """  
         an_extract_annot.py --list_file_annot $annot_file --info_pos $file_rs --out $out --info_file_annot $annot_info
         pdflatex $out
         pdflatex $out
         """
}


bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")

plink_src_ch= Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { plink_src_ch }


fileplotgeno_ch = infors_rs2.join(Channel.from(list_rs).combine(plink_src_ch)).join(Channel.from(list_rs).combine(Channel.fromPath(params.data,checkIfExists:true)))
if(params.covariates)cov_plot_geno="--cov ${params.covariates}"
else  cov_plot_geno=""
if(params.gxe!="")gxe_cov_geno="--gxe ${params.gxe}"
else gxe_cov_geno=""
process PlotByGenotype{
    label 'R'
    memory plink_mem_req
    input :
        set val(rs),file(file_rs), file(bed), file(bim), file(fam), file(data)  from fileplotgeno_ch
    publishDir "${params.output_dir}/$rsnameout", overwrite:true, mode:'copy'
    output :
       set rs, file(outpdf) into report_plot_rs_ch
    script :
        rsnameout=rs.replace(':',"_")
        plinkbase=bed.baseName
        out="plk_"+rsnameout
        outpdf="geno-"+rsnameout+".pdf"
        """
        plink -bfile $plinkbase --extract $file_rs  --recode tab --out $out
        an_plotboxplot.r --ped $out".ped" --data $data --out $outpdf --pheno ${params.pheno} $cov_plot_geno $gxe_cov_geno
        """
}
all_info_rs_ch=report_lz_ch.join(infors_gwas).join(report_info_rs).join(report_plot_rs_ch)
if(params.covariates)covrep="--cov ${params.covariates}"
else  covrep=""
process WriteReportRsFile{
    label 'latex'
    memory plink_mem_req
    input :
       set val(rs), file(locuszoom), file(gwasres),file(annotpdf) ,file(plotgeno) from all_info_rs_ch 
    publishDir "${params.output_dir}/$rsnnameout", overwrite:true, mode:'copy'
    output :
       file(out)
    script :
       rsnnameout=rs.replace(':',"_")
       outtex=rsnnameout+".tex"
       out=rsnnameout+".pdf"
       infoaf=head_freq=="" ? "" :  " --freq_header  $head_freq " 
       """
       an_general_man.py --inp_asso $gwasres --rsname $rs --pheno ${params.pheno} $covrep --out $outtex --chro_header $head_chr --pos_header $head_bp --rs_header $head_rs --pval_header $head_pval --beta_header ${head_beta} $infoaf --locuszoom_plot $locuszoom --geno_plot $plotgeno --inp_asso $gwasres --annot_pdf $annotpdf
       pdflatex $rsnnameout
       pdflatex $rsnnameout
       """ 
}
