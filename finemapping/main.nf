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
 * Description : pipeline to do a finemapping 
 *
 */

def getlistchro(listchro){
 newlistchro=[]
 for(x in listchro.split(',')) {
  splx=x.split("-")
  if(splx.size()==2){
   r1=splx[0].toInteger()
   r2=splx[1].toInteger()
   for(chro in r1..r2){
    newlistchro.add(chro.toString())
   }
  }else if(splx.size()==1){
   newlistchro.add(x)
  }else{
    logger("problem with chro argument "+x+" "+listchro)
    System.exit(0)
  }
 }
 return(newlistchro)
}
//---- General definitions --------------------------------------------------//

import java.nio.file.Paths

// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}
nullfile = [false,"False","false", "FALSE",0,"","0","null",null]
def checkColumnHeader(fname, columns) {
  if (workflow.profile == "awsbatch") return;
  if (fname.toString().contains("s3://")) return;
  if (nullfile.contains(fname)) return;
  new File(fname).withReader { line = it.readLine().tokenize() }
  problem = false;
  columns.each { col ->
    if (! line.contains(col) & col!='') {
      println "The file <$fname> does not contain the column <$col>";
      problem=true;
    }
    if (problem)
      System.exit(2)
  }
}





def helps = [ 'help' : 'help' ]

allowed_params = ["input_dir","input_pat","output","output_dir","data","covariates", "work_dir", "scripts", "max_forks", "cut_maf", "phenotype", "accessKey", "access-key", "secretKey", "secret-key",  "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin", "batch", "batch_col" ,"samplesize", "manifest", "region", "AMI", "queue", "strandreport"]
params_bin=["finemap_bin", "paintor_bin","plink_bin", "caviarbf_bin", "gcta_bin", "gwas_cat_ftp", "list_pheno"]
params_mf=["n_pop","threshold_p", "n_causal_snp", "prob_cred_set"]
params_cojo=["cojo_slct_other", "cojo_top_snps","cojo_slct", "cojo_actual_geno", "threshold_p2", "clump_r2", "size_wind_kb"]
params_filegwas=[ "file_gwas", "head_beta", "head_se", "head_A1", "head_A2", "head_freq", "head_chr", "head_bp", "head_rs", "head_pval", "head_n", "used_pval_z"]
params_paintorcav=["paintor_fileannot", "paintor_listfileannot", "caviarbf_avalue"]
params_memcpu=["gcta_mem_req","plink_mem_req", "plink_cpus_req","other_mem_req","gcta_cpus_req", "fm_cpus_req", "fm_mem_req", "modelsearch_caviarbf_bin","caviar_mem_req"]
param_data=["gwas_cat", "genes_file", "genes_file_ftp"]
param_gccat=["headgc_chr", "headgc_bp", "headgc_bp", "genes_file","genes_file_ftp", "list_chro"]
allowed_params+=params_mf
allowed_params+=params_cojo
allowed_params+=params_filegwas
allowed_params+=params_bin
allowed_params+=params_memcpu
allowed_params+=param_gccat
allowed_params+=params_paintorcav
allowed_params+=param_data



def params_help = new LinkedHashMap(helps)

dummy_dir="${workflow.projectDir}/../qc/input"


params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.genes_file=""
params.genes_file_ftp="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
params.output="finemap"

params.gcta_bin="gcta64"

// paramater
params.n_pop=25000

// params file input
params.head_pval = "P_BOLT_LMM"
params.head_freq = ""
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.head_beta="BETA"
params.head_se="SE"
params.head_A1="ALLELE0"
params.head_A2="ALLELE1"
params.head_n=""
params.used_pval_z=0
params.headgc_chr=""
params.headgc_bp=""
params.gwas_cat = ""

params.prob_cred_set=0.95

params.plink_mem_req="6GB"
params.other_mem_req="20GB"

// gcta parameters
params.gcta_mem_req="15GB"
params.gcta_cpus_req = 1
params.plink_cpus_req=5

params.fm_cpus_req = 5
params.fm_mem_req = "20G"
params.cojo_slct=1
params.cojo_slct_other=""
params.cojo_actual_geno=0
params.big_time='100h'

params.threshold_p=5*10**-8
params.n_causal_snp=1
params.caviarbf_avalue="0.1,0.2,0.4"
params.caviar_mem_req="40GB"
params.paintor_fileannot=""
params.paintor_listfileannot=""
params.threshold_p2=0.5
params.clump_r2=0.5
params.size_wind_kb=100
//params.paintor_annot=""

params.gwas_cat_ftp="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz"
params.list_chro="1-22"
params.list_pheno=""



params.finemap_bin="finemap"
params.caviarbf_bin="caviarbf"
params.modelsearch_caviarbf_bin="model_search"
params.paintor_bin="PAINTOR"
params.plink_bin="plink"

listchro=getlistchro(params.list_chro)
if(params.gwas_cat==""){
println('gwas_cat : gwas catalog option not initialise, will be downloaded')
process GwasCatDl{
    label 'R'
    publishDir "${params.output_dir}/gwascat",  overwrite:true, mode:'copy'
    output :
       file("${out}_all.csv") into gwascat_ch
       file("${out}*")
    script :
      phenol= (params.list_pheno=="") ? "" : "  --pheno '${params.list_pheno}' "
      out="gwascat_format"
      """
      wget -c ${params.gwas_cat_ftp}
      format_gwascat.r --file `basename ${params.gwas_cat_ftp}` $phenol --out $out  --chro ${listchro.join(',')}
      """
}
headgc_chr="chrom"
headgc_bp="chromEnd"
}else{
gwascat_ch=Channel.fromPath(params.gwas_cat, checkIfExists:true)
headgc_chr=params.headgc_chr
headgc_bp=params.headgc_bp
//checkColumnHeader(params.gwas_cat, [headgc_chr,headgc_bp])

}



params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}


bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")

raw_src_ch= Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }


gwas_extract_plk=Channel.create()
plink_subplk=Channel.create()
gwas_plk_clump=Channel.create()
raw_src_ch.separate( gwas_extract_plk, plink_subplk,gwas_plk_clump) { a -> [ a, a,a] }


gwas_file=Channel.fromPath(params.file_gwas,checkIfExists:true)
gwas_file_clump=Channel.fromPath(params.file_gwas,checkIfExists:true)
// plink 
checkColumnHeader(params.file_gwas, [params.head_beta, params.head_se, params.head_A1,params.head_A2, params.head_freq, params.head_chr, params.head_bp, params.head_rs, params.head_pval, params.head_n])
process clump_data{
 cpus params.plink_cpus_req
 input :
     set file(bed),file(bim),file(fam) from gwas_plk_clump
     file(gwasfile) from gwas_file_clump
 publishDir "${params.output_dir}/clump/", overwrite:true, mode:'copy'
 output :
    file("${output}.clumped") into file_clump 
 script :
   output="clump_output"
   plkfile=bed.baseName
   """
   plink_format.py --inp_asso $gwasfile --chro_header  ${params.head_chr} --bp_header ${params.head_bp} --a1_header ${params.head_A1} --a2_header ${params.head_A2}  --pval_header ${params.head_pval} --beta_header ${params.head_beta}  --out $output --rs_header ${params.head_rs} --se_header  ${params.head_se}
 ${params.plink_bin} -bfile $plkfile  -out $output --keep-allele-order --threads ${params.plink_cpus_req}   --clump $output --clump-p1 ${params.threshold_p} --clump-p2 ${params.threshold_p2} --clump-r2 ${params.clump_r2} --clump-kb ${params.size_wind_kb}
  """
}
process extract_sigpos{
  input :
   file(clump) from file_clump
  output :
     stdout into (postonalyse, postonalyse2)
  script:
  """
   sed '1d' $clump| awk '{if(\$1!="")print \$1"_"\$4}' 
  """

}
pos_toanalyse_ch = Channel.create()
check = Channel.create()

//  chrolist=chrolist.flatMap { list_str -> list_str.split() }
pos_tonalyse_ch=postonalyse.flatMap { list_str -> list_str.split() }


process ExtractPositionGwas{
  memory params.other_mem_req
  input :
     file(filegwas) from gwas_file
     set file(bed),file(bim),file(fam) from gwas_extract_plk
  each  pos from pos_tonalyse_ch
  output :
    set val(pos),file("${out}.gcta") into gcta_gwas
    set val(pos),file("${out}_finemap.z") into  (finemap_gwas_cond, finemap_gwas_sss)
    set val(pos),file("${out}_caviar.z") into caviarbf_gwas
    set val(pos),file("${out}.paintor") into paintor_gwas
    set val(pos),file("${out}.range") into range_plink
    set val(pos), file("${out}.all") into data_i
    set val(pos), file("${out}.pos") into paintor_gwas_annot
  publishDir "${params.output_dir}/$pos/file_format/", overwrite:true, mode:'copy'
  script :
    bp=pos.split('_')[1]
    chro=pos.split('_')[0]
    freq= (params.head_freq=="") ? "":" --freq_header ${params.head_freq} "
    nheader= (params.head_n=="") ? "":" --n_header ${params.head_n}"
    nvalue= (params.n_pop=="") ? "":" --n ${params.n_pop}"
    out=pos
    bfile=bed.baseName
    around=params.size_wind_kb*1000 
    """
    end_seq=`expr $bp + $around`
    begin_seq=`expr $bp - $around`
    fine_extract_sig.py --inp_resgwas $filegwas --chro ${chro} --begin \$begin_seq  --end \$end_seq --chro_header ${params.head_chr} --pos_header ${params.head_bp} --beta_header ${params.head_beta} --se_header ${params.head_se} --a1_header ${params.head_A1} --a2_header ${params.head_A2} $freq  --bfile $bfile --rs_header ${params.head_rs} --out_head $out --p_header ${params.head_pval}  $nvalue --min_pval ${params.threshold_p} $nheader --z_pval ${params.used_pval_z}
    """
}

range_plink_ch=range_plink.combine(plink_subplk)

process SubPlink{
  input :
     set val(pos),file(range), file(bed),file(bim),file(fam) from range_plink_ch
  output :
     set val(pos),file("${out}.bed"),file("${out}.bim"),file("${out}.fam") into (subplink_ld, subplink_gcta)
  script : 
     plk=bed.baseName
     out=plk+'_sub'
     """
     ${params.plink_bin} -bfile $plk  --keep-allele-order --extract  range  $range --make-bed -out  $out
     """
}

process ComputedLd{
   memory params.other_mem_req
   input : 
      set val(pos),file(bed),file(bim),file(fam) from subplink_ld
  output :
       set val(pos),file("$outld") into (ld_fmcond, ld_fmsss,ld_caviarbf, ld_paintor)
   script :
    outld=pos+".ld"
    plk=bed.baseName
    """
    ${params.plink_bin} --r2 square0 yes-really -bfile $plk -out "tmp"
    sed 's/\\t/ /g' tmp.ld | sed 's/nan/0/g' > $outld
    """
}

ld_fmcond_group=ld_fmcond.join(finemap_gwas_cond)

process ComputedFineMapCond{
  label 'finemapping'
  cpus params.fm_cpus_req
  memory params.fm_mem_req
  input :
    set val(pos),file(ld),file(filez) from ld_fmcond_group
  publishDir "${params.output_dir}/$pos/fm_cond", overwrite:true, mode:'copy'
  output :
    set val(pos), file("${out}.snp") into res_fmcond
    set file("${out}.config"), file("${out}.cred"), file("${out}.log_cond")
  script:
  fileconfig="config"
  out=pos+"_cond" 
  """ 
  echo "z;ld;snp;config;cred;log;n_samples" > $fileconfig
  echo "$filez;$ld;${out}.snp;${out}.config;${out}.cred;${out}.log;${params.n_pop}" >> $fileconfig
  ${params.finemap_bin} --cond --in-files $fileconfig   --log --cond-pvalue ${params.threshold_p}  --n-causal-snps ${params.n_causal_snp}  --prob-cred-set ${params.prob_cred_set}
  """
}

ld_fmss_group=ld_fmsss.join(finemap_gwas_sss)

process ComputedFineMapSSS{
  label 'finemapping'
  memory params.fm_mem_req
  cpus params.fm_cpus_req
  input :
     set val(pos),file(ld),file(filez) from ld_fmss_group
  publishDir "${params.output_dir}/$pos/fm_sss", overwrite:true, mode:'copy'
  output :
    set val(pos),file("${out}.snp") into res_fmsss
    set file("${out}.config"), file("${out}.cred${params.n_causal_snp}"), file("${out}.log_sss")
  script:
  fileconfig="config"
  out=pos+"_sss"
  """
  echo "z;ld;snp;config;cred;log;n_samples" > $fileconfig
  echo "$filez;$ld;${out}.snp;${out}.config;${out}.cred;${out}.log;${params.n_pop}" >> $fileconfig
  ${params.finemap_bin} --sss --in-files $fileconfig  --n-threads ${params.fm_cpus_req}  --log --n-causal-snps ${params.n_causal_snp} --prob-cred-set ${params.prob_cred_set}
  """
}

ld_caviarbf_group=ld_caviarbf.join(caviarbf_gwas)
process ComputedCaviarBF{
  memory params.caviar_mem_req
  label 'finemapping'
  input :
     set val(pos),file(ld),file(filez) from ld_caviarbf_group
  publishDir "${params.output_dir}/$pos/caviarbf", overwrite:true, mode:'copy'
  output :
   set val(pos), file("${output}.marginal") into res_caviarbf
   set file("$output"), file("${output}.statistics")
  script :
   output=pos+"_caviarbf"
   """
   ${params.caviarbf_bin} -z ${filez} -r $ld  -t 0 -a ${params.caviarbf_avalue} -c ${params.n_causal_snp} -o ${output} -n ${params.n_pop}
   nb=`cat ${filez}|wc -l `
   ${params.modelsearch_caviarbf_bin} -i $output -p 0 -o $output -m \$nb 2> ${output}_modelsearch.log
   """
}

//NCausalSnp=Channel.from(1..params.n_causal_snp)
if(params.paintor_bin!="0" & params.paintor_bin!=0){
baliseannotpaint=0
if(params.paintor_fileannot!=""){
paintor_fileannot=Channel.fromPath(params.paintor_fileannot)
paintor_fileannotplot=Channel.fromPath(params.paintor_fileannot)
baliseannotpaint=1
}else{
 if(params.paintor_listfileannot!=""){
  baliseannotpaint=1
  paintor_gwas_annot2=paintor_gwas_annot.combine(Channel.fromPath(params.paintor_listfileannot))
  process paintor_selectannot{
   input :
    set val(pos),file(list_loc), file(listinfo) from paintor_gwas_annot2
   publishDir "${params.output_dir}/$pos/paintor/annot", overwrite:true, mode:'copy'
   output :
    set val(pos),file(out) into (paintor_fileannot, paintor_fileannotplot, paintor_fileannot2)
   script :
   outtmp="tmp.res"
   out="annotationinfo"
   """
   head -1 $list_loc > $outtmp
   sed '1d' $list_loc |awk '{print "chr"\$0}' >> $outtmp
   annotate_locus_paint.py --input $listinfo  --locus $outtmp --out $out --chr chromosome --pos position
   """
  }
  process paintor_extractannotname{
    input :
       set val(pos),file(fileannot) from paintor_fileannot2
    output :
       set val(pos),stdout into annotname
    """
    head -1 $fileannot | sed 's/ /,/g' 
    """ 
  }
} else{
println 'no file annot for paintor'
postonalyse2.into{postonanalyse_tmp1; postonanalyse_tmp2;postonanalyse_tmp3}

pos_tonalyse_ch_1=postonanalyse_tmp1.flatMap { list_str -> list_str.split() }
paintor_fileannot=pos_tonalyse_ch_1.combine(Channel.fromPath("${dummy_dir}/0"))

pos_tonalyse_ch_2=postonanalyse_tmp2.flatMap { list_str -> list_str.split() }
paintor_fileannotplot=pos_tonalyse_ch_2.combine(Channel.fromPath("${dummy_dir}/0"))

pos_tonalyse_ch_3=postonanalyse_tmp3.flatMap { list_str -> list_str.split() }
annotname=pos_tonalyse_ch_3.combine(Channel.value("N"))
}
}


ld_paintor_group=ld_paintor.join(paintor_gwas).join(paintor_fileannot).join(annotname) //.combine(annotname).combine(paintor_fileannot)
process ComputedPaintor{
   label 'finemapping'
   memory params.fm_mem_req
   input :
    set val(pos),file(ld),file(filez), file(fileannot), val(annot_name) from ld_paintor_group
  publishDir "${params.output_dir}/$pos/paintor/", overwrite:true, mode:'copy'
  output :
      set val(pos),file("${output}.results") into res_paintor_ch
      set val(pos),file(FileInfo) into infores_paintor_ch
      file("${output}*")
      set val(pos), file(BayesFactor) into res_paintor_ch_bf
  script :
    ncausal=params.n_causal_snp
    output=pos+"_paintor_$ncausal" 
    DirPaintor=output
    annot=(baliseannotpaint==0) ? "" : " -Gname ${output}_an  -annotations ${annot_name}"
    BayesFactor=output+".BayesFactor"
    FileInfo=output+".info"
    Info="$ncausal;${output}.results;$BayesFactor"
    """
    echo "$Info" > $FileInfo
    echo $output > input.files
    cp $filez $output
    cp $ld $output".ld"
    if [ $fileannot == "0" ]
    then
    paint_annotation.py $fileannot $output  $output".annotations"
    else 
    cp $fileannot $output".annotations"
    fi
    ${params.paintor_bin} -input input.files -in ./ -out ./ -Zhead Z -LDname ld -enumerate $ncausal -num_samples  ${params.n_pop} -Lname $BayesFactor $annot
    """
}
}else{
 ld_paintor_group=ld_paintor.join(paintor_gwas)
 process ComputedPaintor_null{
  input :
   set val(pos), file(ld),file(filez) from ld_paintor_group
  output :
     set val(pos),file("${output}.results") into res_paintor_ch
     set val(pos),file(FileInfo) into infores_paintor_ch
  script :
     output=pos+"_paintor_null"
     FileInfo=output+".info"
     BayesFactor=output+".BayesFactor"
     """
     touch $output".results"
     touch $BayesFactor
     touch $FileInfo
     """
 }
paintor_fileannot=file("${dummy_dir}/0")
paintor_fileannotplot=file("${dummy_dir}/0")
annotname=Channel.from("N")

}
//res_paintor_ch=res_paintor.collect()
//infores_paintor_ch=infores_paintor.collect()

gcta_gwas_join=gcta_gwas.join(subplink_gcta)
process ComputedCojo{
   label 'gcta'
   memory params.gcta_mem_req
   cpus params.gcta_cpus_req
   input :
     set val(pos),file(filez), file(bed),file(bim),file(fam) from gcta_gwas_join
   publishDir "${params.output_dir}/$pos/cojo_gcta", overwrite:true, mode:'copy'
   output :
     set val(pos), file("${output}.jma.cojo")  into res_cojo
     set file("${output}.cma.cojo"), file("${output}.ldr.cojo"), file("${output}.log")
   script :
    output=pos+"_cojo"
    plk=bed.baseName
    """ 
    ${params.gcta_bin} --bfile $plk  --cojo-slct --cojo-file $filez --out $output  --cojo-p ${params.threshold_p} --thread-num ${params.gcta_cpus_req}  --diff-freq 0.49
    if [ ! -f $output".jma.cojo" ]
    then 
    echo -e  "Chr	SNP	bp	refA	freq	b	se	p	n	freq_geno	bJ	bJ_se	pJ	LD_r"	> $output".jma.cojo"
    touch $output".cma.cojo" $output".ldr.cojo"
    fi
    """

}
if(params.genes_file==""){
process GetGenesInfo{
   label 'R'
   output :
      file(out) into genes_file_ch
   publishDir "${params.output_dir}/data/", overwrite:true, mode:'copy'
   script :
     out="gencode.v19.genes"
     """
     wget -c ${params.genes_file_ftp}
     zcat `basename ${params.genes_file_ftp}` > file_genes
     change_genes.r file_genes
     """
}
}else{
genes_file_ch=Channel.fromPath(params.genes_file)
}

mergeall=res_paintor_ch.join(infores_paintor_ch).join(res_paintor_ch_bf).join(paintor_fileannotplot).join(res_cojo).join(res_caviarbf).join(res_fmsss).join(res_fmcond).join(data_i).combine(genes_file_ch).combine(gwascat_ch)
process MergeResult{
    label 'R'
    memory params.other_mem_req
    input :
      set val(pos), file(paintor),file(infopaintor), file(paintor_bf), file(pfileannot), file(cojo), file(caviarbf), file(fmsss), file(fmcond), file(datai), file(genes), file(gwascat) from mergeall
   publishDir "${params.output_dir}/$pos/", overwrite:true, mode:'copy'
    output :
       set file("${out}.pdf"), file("${out}.all.out"), file("${out}.all.out")
    script :
      out=params.output
      //infopaint=infopaintor.join(" ")
      pfileannot= (baliseannotpaint=="0" | params.paintor_bin==0) ? "":" --paintor_fileannot $pfileannot "
      paintor = (params.paintor_bin=="0") ? "" : "--listpaintor  infopaintor"
      
      """
       cat $infopaintor > infopaintor
       echo "sss $fmsss" > infofinemap 
       echo "cond $fmcond" >> infofinemap 
       merge_finemapping_v2.r --out $out --listpaintor  infopaintor  --cojo  $cojo --datai  $datai --caviarbf $caviarbf --list_genes $genes  --gwascat $gwascat --headbp_gc ${headgc_bp} --headchr_gc ${headgc_chr}  --listfinemap infofinemap  $pfileannot
      """

}

