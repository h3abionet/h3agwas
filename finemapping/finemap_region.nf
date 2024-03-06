#!/usr/bin/env nextflow

/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2022
 *
 *
 * Description : pipeline to do a finemapping 
 *
 */

nextflow.enable.dsl = 1
if (!workflow.resume) {
    def dir = new File(params.output_dir)
    if (dir.exists() && dir.directory && (!(dir.list() as List).empty)) {
       println "\n\n============================================"
       println "Unless you are doing a -resume, the output directory should be empty"
       println "We do not want to overwrite something valuable in "+params.output_dir
       println "Either clean your output directory or check if you meant to do a -resume"
       System.exit(-1)
    }
}


def strmem(val){
 return val as nextflow.util.MemoryUnit
}


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
    System.exit(-1)
  }
 }
 return(newlistchro)
}

/*definition*/
def errormess(message,exitn=0){
    if(message=="")return(0)
    println(message)
    System.exit(exitn)
}
def checkparams(param, namesparam, type, min=null, max=null, possibleval=null, notpossibleval=null) {
  messageerror=""
  if(param==null){
    messageerror+="error :--"+namesparam+" is null "
  } else {
    if(!(param.getClass() in type)){
      messageerror+="error :--"+namesparam+" must be a "+ type
     if(params.getClass()==Boolean)messageerror+=", but no parameters given"
     else messageerror+=" but type is "+param.getClass()+" value "+ param
   }else{
   if(min && param<min)messageerror+="\nerror : --"+namesparam+" < min value :"+param +" < "+min
   if(max && param>max)messageerror+="\nerror : --"+namesparam +"> maxvalue :" + param+" > "+max
   if(possibleval && !(param in possibleval))messageerror+="\nerro : --"+namesparam +" must be one the value :"+possibleval.join(',')
   }
   }
    errormess(messageerror,-1)
}


def checkmultiparam(params, listparams, type, min=null, max=null, possibleval=null, notpossibleval=null){
 messageerror=""
 for(param in listparams){
   if(params.containsKey(param)){
     checkparams(params[param], param, type, min=min, max=max, possibleval=possibleval, notpossibleval=notpossibleval)
   }else{
     messageerror+="param :"+param+" unknown (check your command line)\n"
   }
 }
 errormess(messageerror, -1)
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
  if (fname.toString().contains("az://")) return;
  if (nullfile.contains(fname)) return;
  new File(fname).withReader { line = it.readLine().tokenize() }
  problem = false;
  columns.each { col ->
    if (! line.contains(col) & col!='') {
      println "The file <$fname> does not contain the column <$col>";
      problem=true;
    }
    if (problem)
      System.exit(-1)
  }
}





def helps = [ 'help' : 'help' ]


allowed_params_input = ["input_dir","input_pat","output","output_dir","plink_mem_req", "work_dir", "scripts",  "accessKey", "access-key", "secretKey", "secret-key", "region",  "big_time",  "rs_list", 'list_phenogc', 'cojo_slct_other', "paintor_fileannot", "paintor_listfileannot", "caviarbf_avalue", "gwas_cat", "genes_file", "genes_file_ftp", "list_phenogc", "file_phenogc", "headgc_chr", "headgc_bp", "headgc_bp", "genes_file","genes_file_ftp", "list_chro", 'file_pheno', 'modelsearch_caviarbf_bin', "AMI", "instanceType", "instance-type", "bootStorageSize","maxInstances", "max-instances", "sharedStorageMount", "shared-storage-mount",'queue',"ftp_vcf", 'file_gwas',  "list_vcf"]
allowed_params=allowed_params_input
allowed_params_bin=["finemap_bin", "paintor_bin","plink_bin", "caviarbf_bin", "gcta_bin", "gwas_cat_ftp"]
allowed_params+=allowed_params_bin
allowed_params_cores=["plink_cpus_req", "gcta_cpus_req", "fm_cpus_req",'max_plink_cores', "other_cpus_req"]
allowed_params+=allowed_params_cores
allowed_params_intother=["max_forks", "n_pop", "n_causal_snp", 'cojo_slct', "size_wind_kb", "begin_seq", "end_seq"]
allowed_params+=allowed_params_intother
allowed_params_bolother=["used_pval_z", "ftp1000genome"]
allowed_params+=allowed_params_bolother
allowed_params_float=["cut_maf", "threshold_p", "prob_cred_set"]
allowed_params+=allowed_params_float
allowed_params_memory=["gcta_mem_req" , "plink_mem_req", "other_mem_req", "fm_mem_req","caviar_mem_req", "paintor_mem_req", "boot-storage-size"]
allowed_params+=allowed_params_memory
params_filegwas=[ "head_beta", "head_se", "head_A1", "head_A2", "head_freq", "head_chr", "head_bp", "head_rs", "head_pval", "head_n"]
allowed_params+=params_filegwas
params_strorint=["chro"]
allowed_params+=params_strorint




def params_help = new LinkedHashMap(helps)
filescript=file(workflow.scriptFile)
projectdir="${filescript.getParent()}"
dummy_dir="${projectdir}/../qc/input"


params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.genes_file=""
params.genes_file_ftp="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
params.output="finemap"
params.file_phenogc=""


// paramater
params.n_pop=0

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
params.ftp1000genome=1
params.list_vcf=""
params.cut_maf=0.01

params.prob_cred_set=0.95


params.other_mem_req="20GB"

// gcta parameters
params.gcta_mem_req="15GB"
params.gcta_cpus_req = 1

params.fm_cpus_req = 5
params.cojo_slct=1
params.cojo_slct_other=""
params.big_time='100h'

params.threshold_p=5*10**-8
params.n_causal_snp=3
params.caviarbf_avalue="0.1,0.2,0.4"

params.caviar_mem_req="40GB"
params.paintor_mem_req="20GB"
params.fm_mem_req = "20G"
params.plink_mem_req="6GB"

params.paintor_fileannot=""
params.paintor_listfileannot=""
params.other_cpus_req=10
//params.paintor_annot=""

params.gwas_cat_ftp="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz"
params.list_phenogc=""
params.ftp_vcf="ftp://ftp.1000genomes.ebi.ac.uk:21/vol1/ftp/release/20130502/"



params.finemap_bin="finemap"
params.caviarbf_bin="caviarbf"
params.modelsearch_caviarbf_bin="model_search"
params.paintor_bin="PAINTOR"
params.plink_bin="plink"


params.chro=""
params.begin_seq=""
params.end_seq=""

checkmultiparam(params,allowed_params_input, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_memory, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_bin, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_cores, java.lang.Integer, min=1, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_intother, java.lang.Integer, min=0, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_bolother, java.lang.Integer, min=0, max=null, possibleval=[0,1], notpossibleval=null)
checkmultiparam(params,allowed_params_float, [java.lang.Double,java.lang.Float, java.lang.Integer, java.math.BigDecimal], min=0, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,params_filegwas, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,params_strorint, [java.lang.String,java.lang.Integer], min=null, max=null, possibleval=null, notpossibleval=null)


if(params.begin_seq > params.end_seq){
error('begin_seq > end_seq')
}
if(params.input_dir=="" | params.input_pat==""){
 
 if(params.ftp1000genome==1){
  listchro_ch=channel.from(params.chro)
  process Dl1000G{
   label 'py3utils'
   cpus params.other_cpus_req
   input :
      val(chro) from listchro_ch
   output :
       tuple val(chro), file("${file1000G}") into file1000G
   script :
      //ftp://ftp.1000genomes.ebi.ac.uk:21/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      file1000G= (chro=='X') ? "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz" : "ALL.chr${chro}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
      file1000G= (chro=='Y') ? "ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz" : "$file1000G"
      """
      echo  $chro $params.begin_seq $params.end_seq | awk '{print \$1"\t"\$2"\t"\$3"\t"\$3":"\$2}'> pos.bed
      bcftools view --threads ${params.other_cpus_req} -R pos.bed ${params.ftp_vcf}/$file1000G|bgzip -c > $file1000G
      """
  }
 }else{
  if(params.list_vcf!=''){
    listvcf_ch=Channel.fromPath(file(params.list_vcf, checkIfExists:true).readLines(), checkIfExists:true)
    process findchro{
      label 'py3utils'
      input :
        path(vcf) from listvcf_ch
      output :
        tuple env(chro), path(vcf), path("${vcf}.tbi") into (filevcf_chro_ch1, filevcf_chro_ch2)
      script :
        """
        chro=`zcat $vcf|grep -v "#"|head -1|awk '{print \$1}'`
        tabix -p vcf $vcf
        """
     }
   filevcf_chro_charray=filevcf_chro_ch1.join(Channel.from(params.chro))
   process extract_chrovcf{
     label 'py3utils'
     input :
        tuple val(chro),path(vcf), path(vcfindex) from filevcf_chro_charray
      output :
        tuple val(chro),  path(vcfout) into file1000G
      script :
        vcfout="vcf_array_"+chro+".bcf"
        """
        awk -v chro=$chro '{if(\$1==chro)print \$1"\\t"\$2-1"\\t"\$2"\\t"\$1":"\$2}' $pos_geno > pos.bed
        bcftools view --threads ${params.other_cpus_req} -R pos.bed $vcf|bgzip -c > $vcfout
        """
   }

  }else{
   println "error not params.list_vcf"
  }
 }

 process transfvcfInBed1000G{
   cpus params.other_cpus_req
   input :
     tuple val(chro), file(vcf1000) from file1000G
   output :
     tuple val(chro), path("${out}.bim"), path("${out}.fam"), path("${out}.bed") into plk_chro
   script :
     out="chrtmp_"+chro
     """
     plink --vcf $vcf1000  --keep-allele-order  --make-bed -out $out --threads ${params.other_cpus_req}
     """
 }
//
 process cleanPlinkFile{
    cpus params.other_cpus_req
    input :
     tuple val(chro), path(bim), path(fam), path(bed) from plk_chro
    output :
     tuple path("${out}.bim"), path("${out}.fam"), path("${out}.bed") into plk_chro_cl
    script :
    plk=bim.baseName
    out="chr_"+chro+"_clean"
      """
      cp $fam $plk"_tmp.fam"
      cp $bed $plk"_tmp.bed"
      awk '{if(\$2=="."){\$2=\$1":"\$4};print \$0}' "${bim}" > $plk"_tmp.bim"
      awk '{if(length(\$5)==1 && length(\$6)==1)print \$2}' $plk"_tmp.bim" > ${bim}.wellpos.pos
      awk '{print \$2}' $plk"_tmp.bim" | sort | uniq -d > duplicated_snps.snplist
      plink -bfile $plk"_tmp"  --keep-allele-order --extract  ${bim}.wellpos.pos --make-bed -out $out --threads ${params.other_cpus_req} --exclude duplicated_snps.snplist 
     """
 }
 plk_chro_flt=plk_chro_cl.collect()

 process mergePlinkFile{
   cpus params.other_cpus_req
   input :
     file(allfile) from plk_chro_flt
   output :
     tuple file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into raw_src_ch
   script :
     allfile2=allfile.toList().collect{it.toString().replaceFirst(~/\.[^\.]+$/, '')}
     allfile2=allfile2.unique()
     firstbed=allfile2[0]
     allfile2.remove(0)
     nbfile = allfile2.size()
     out=params.output+"_befsex"
     """
     if [ $nbfile == "0" ]
     then
     cp ${firstbed}.fam ${out}.fam
     cp ${firstbed}.bed ${out}.bed
     cp ${firstbed}.bim ${out}.bim
     else
     echo "${allfile2.join('\n')}" > listbedfile
     plink --bfile $firstbed  --keep-allele-order --threads ${params.other_cpus_req} --merge-list listbedfile --make-bed --out $out -maf ${params.cut_maf}
     fi

     """
 }
 
}else{

bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")

raw_src_ch= Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }
}




if(params.gwas_cat==""){
println('gwas_cat : gwas catalog option not initialise, will be downloaded')
if(params.file_phenogc=="")phenogc_ch=channel.fromPath("${dummy_dir}/01")
else phenogc_ch=channel.fromPath(params.file_phenogc,checkIfExists:true)


process GwasCatDl{
    label 'R'
    publishDir "${params.output_dir}/gwascat",   mode:'copy'
    input :
        file(phenogc) from  phenogc_ch
    output :
       file("${out}_all.csv") into gwascat_ch
       file("${out}*")
    script :
      phenol= (params.list_phenogc=="") ? "" : "  --pheno '${params.list_phenogc}' "
      phenofile= (params.file_phenogc=="") ? "" : "  --file_pheno  $phenogc "
      out="gwascat_format"
      """
      wget -c ${params.gwas_cat_ftp} --no-check-certificate
      format_gwascat.r --file `basename ${params.gwas_cat_ftp}` $phenol --out $out  --chro ${params.chro} $phenofile
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
if(params.chro=="" | params.begin_seq=="" | params.end_seq==""){
error('chro, begin_seq or end_seq not initialise')
}



params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}



gwas_extract_plk=Channel.create()
plink_subplk=Channel.create()
raw_src_ch.separate( gwas_extract_plk, plink_subplk) { a -> [ a, a] }


gwas_file=Channel.fromPath(params.file_gwas,checkIfExists:true)
// plink 
checkColumnHeader(params.file_gwas, [params.head_beta, params.head_se, params.head_A1,params.head_A2, params.head_freq, params.head_chr, params.head_bp, params.head_rs, params.head_pval, params.head_n])
process ExtractPositionGwas{
  memory params.plink_mem_req
  cpus params.max_plink_cores
  input :
     file(filegwas) from gwas_file
     set file(bed),file(bim),file(fam) from gwas_extract_plk
  output :
    file("${out}.gcta") into gcta_gwas
    file("${out}_finemap.z") into  (finemap_gwas_cond, finemap_gwas_sss)
    file("${out}_caviar.z") into caviarbf_gwas
    file("${out}.paintor") into paintor_gwas
    file("${out}.range") into range_plink
    file("${out}.all") into data_i
    file("${out}.pos") into paintor_gwas_annot
    env n into nval_ch, nval_ch_fmss, nval_ch_fm, nval_ch_caviar, nval_ch_paintor
  publishDir "${params.output_dir}/file_format/",  mode:'copy'
  script :
    freq= (params.head_freq=="") ? "":" --freq_header ${params.head_freq} "
    nheader= (params.head_n=="") ? "":" --n_header ${params.head_n}"
    nvalue= (params.n_pop<=1) ? "":" --n ${params.n_pop}"
    out=params.chro+"_"+params.begin_seq+"_"+params.end_seq
    bfile=bed.baseName
    plink_mem_req_max=params.plink_mem_req.replace('GB','000').replace('KB','').replace(' ','').replace('MB','').replace('Mb','')
    """
    fine_extract_sig.py --inp_resgwas $filegwas --chro ${params.chro} --begin ${params.begin_seq}  --end ${params.end_seq} --chro_header ${params.head_chr} --pos_header ${params.head_bp} --beta_header ${params.head_beta} --se_header ${params.head_se} --a1_header ${params.head_A1} --a2_header ${params.head_A2} $freq  --bfile $bfile --rs_header ${params.head_rs} --out_head $out --p_header ${params.head_pval}  $nvalue --min_pval ${params.threshold_p} $nheader --z_pval ${params.used_pval_z} --maf ${params.cut_maf} --threads ${params.max_plink_cores} 
    n=`cat n.out`
    """
}


process SubPlink{
  memory params.plink_mem_req
  cpus params.max_plink_cores
  input :
     set file(bed),file(bim),file(fam) from plink_subplk
     file(range) from range_plink
  output :
     set file("${out}.bed"),file("${out}.bim"),file("${out}.fam") into (subplink_ld, subplink_gcta)
  script : 
     plk=bed.baseName
     out=plk+'_sub'
    plink_mem_req_max=params.plink_mem_req.replace('GB','000').replace('KB','').replace(' ','').replace('MB','').replace('Mb','')
     """
     ${params.plink_bin} -bfile $plk  --keep-allele-order --extract  range  $range --make-bed -out  $out 
     """
}

process ComputedLd{
   memory params.plink_mem_req
   cpus params.max_plink_cores
   input : 
      set file(bed),file(bim),file(fam) from subplink_ld
  output :
       file("$outld") into (ld_fmcond, ld_fmsss,ld_caviarbf, ld_paintor)
   script :
    outld=params.chro+"_"+params.begin_seq+"_"+params.end_seq+".ld"
    plk=bed.baseName
  plink_mem_req_max=params.plink_mem_req.replace('GB','000').replace('KB','').replace(' ','').replace('MB','').replace('Mb','')

    """
    ${params.plink_bin} --r2 square0 yes-really -bfile $plk -out "tmp" --threads ${params.max_plink_cores}  --memory  $plink_mem_req_max
    sed 's/\\t/ /g' tmp.ld | sed 's/nan/0/g' > $outld
    """
}


process ComputedFineMapCond{
  label 'finemapping'
  cpus 2
  memory params.fm_mem_req
  input :
    file(ld) from ld_fmcond 
    file(filez) from finemap_gwas_cond
    val(n) from nval_ch_fm
  publishDir "${params.output_dir}/fm_cond",  mode:'copy'
  output :
    tuple path("${out}.snp"), path("${out}.cred") into res_fmcond
    set file("${out}.config"), file("${out}.cred"), file("${out}.log_cond")
  script:
  fileconfig="config"
  out=params.chro+"_"+params.begin_seq+"_"+params.end_seq+"_cond" 
  """ 
  echo "z;ld;snp;config;cred;log;n_samples" > $fileconfig
  echo "$filez;$ld;${out}.snp;${out}.config;${out}.cred;${out}.log;$n" >> $fileconfig
  ${params.finemap_bin} --cond --in-files $fileconfig   --log --cond-pvalue ${params.threshold_p}  --n-causal-snps ${params.n_causal_snp}  --prob-cred-set ${params.prob_cred_set} 
  """
}

process ComputedFineMapSSS{
  label 'finemapping'
  memory params.fm_mem_req
  cpus params.fm_cpus_req
  input :
    file(ld) from ld_fmsss
    file(filez) from finemap_gwas_sss
    val(n) from nval_ch_fmss
  publishDir "${params.output_dir}/fm_sss",  mode:'copy'
  output :
    tuple path("${out}.snp"), path("${out}.cred${params.n_causal_snp}") into res_fmsss
    set file("${out}.config"), file("${out}.cred${params.n_causal_snp}"), file("${out}.log_sss")
  script:
  fileconfig="config"
  out=params.chro+"_"+params.begin_seq+"_"+params.end_seq+"_sss"
  """
  echo "z;ld;snp;config;cred;log;n_samples" > $fileconfig
  echo "$filez;$ld;${out}.snp;${out}.config;${out}.cred;${out}.log;$n" >> $fileconfig
  ${params.finemap_bin} --sss --in-files $fileconfig  --n-threads ${params.fm_cpus_req}  --log --n-causal-snps ${params.n_causal_snp} --prob-cred-set ${params.prob_cred_set}
  """
}

process ComputedCaviarBF{
  memory params.caviar_mem_req
  label 'finemapping'
  input :
    file(filez) from caviarbf_gwas
    file(ld) from ld_caviarbf
    val(n) from nval_ch_caviar
  publishDir "${params.output_dir}/caviarbf",  mode:'copy'
  output :
   file("${output}.marginal") into res_caviarbf
   set file("$output"), file("${output}.statistics")
  script :
   output=params.chro+"_"+params.begin_seq+"_"+params.end_seq+"_caviarbf"
   """
   ${params.caviarbf_bin} -z ${filez} -r $ld  -t 0 -a ${params.caviarbf_avalue} -c ${params.n_causal_snp} -o ${output} -n $n
   nb=`cat ${filez}|wc -l `
   ${params.modelsearch_caviarbf_bin} -i $output -p 0 -o $output -m \$nb 2> ${output}_modelsearch.log
   """
}

NCausalSnp=Channel.from(1..params.n_causal_snp)
baliseannotpaint=0
if(params.paintor_fileannot!=""){
paintor_fileannot=Channel.fromPath(params.paintor_fileannot)
paintor_fileannotplot=Channel.fromPath(params.paintor_fileannot)
baliseannotpaint=1
}else{
 if(params.paintor_listfileannot!=""){
  baliseannotpaint=1
  paintor_listfileannot=Channel.fromPath(params.paintor_listfileannot)
  process paintor_selectannot{
   input :
    file(listinfo) from paintor_listfileannot
    file(list_loc) from paintor_gwas_annot
   publishDir "${params.output_dir}/paintor/annot",  mode:'copy'
   output :
    file(out) into (paintor_fileannot, paintor_fileannotplot, paintor_fileannot2)
   script :
   outtmp="tmp.res"
   out="annotationinfo"
   """
   head -1 $list_loc > $outtmp
   sed '1d' $list_loc |awk '{print "chr"\$0}' >> $outtmp
   annotate_locus_paint.py --input $listinfo  --locus $outtmp --out $out --chr chromosome --pos position
   """
  }
  paintor_listfileannot2=Channel.fromPath(params.paintor_listfileannot)
  process paintor_extractannotname{
    input :
       file(fileannot) from paintor_fileannot2
    output :
       stdout into annotname
    """
    head -1 $fileannot | sed 's/ /,/g' 
    """ 
  }
} else{
paintor_fileannot=file("${dummy_dir}/0")
paintor_fileannotplot=file("${dummy_dir}/0")
annotname=Channel.from("N")
}
}
process ComputedPaintor{
   label 'finemapping'
   memory params.paintor_mem_req
   input :
    file(filez) from paintor_gwas
    file(ld) from ld_paintor
    file(fileannot) from paintor_fileannot
    val(annot_name) from annotname
    val(n) from nval_ch_paintor
  each ncausal from NCausalSnp
  publishDir "${params.output_dir}/paintor/",  mode:'copy'
  output :
      set file("${output}.results"), file("$BayesFactor") into res_paintor
      file(FileInfo) into infores_paintor
      file("${output}*")
  script :
    output=params.chro+"_"+params.begin_seq+"_"+params.end_seq+"_paintor_$ncausal" 
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
    ${params.paintor_bin} -input input.files -in ./ -out ./ -Zhead Z -LDname ld -enumerate $ncausal -num_samples  $n -Lname $BayesFactor $annot
    """
}
res_paintor_ch=res_paintor.collect()
infores_paintor_ch=infores_paintor.collect()


process ComputedCojo{
   label 'gcta'
   memory params.gcta_mem_req
   cpus params.gcta_cpus_req
   input :
     set  file(bed),file(bim),file(fam) from subplink_gcta
     file(filez) from gcta_gwas
   publishDir "${params.output_dir}/cojo_gcta",  mode:'copy'
   output :
     file("${output}.jma.cojo")  into res_cojo
     set file("${output}.cma.cojo"), file("${output}.ldr.cojo"), file("${output}.log")
   script :
    output=params.chro+"_"+params.begin_seq+"_"+params.end_seq+"_cojo"
    plk=bed.baseName
    """ 
    ${params.gcta_bin} --bfile $plk  --cojo-slct --cojo-file $filez --out $output  --cojo-p ${params.threshold_p} --thread-num ${params.gcta_cpus_req}  --diff-freq 0.49
    """

}
if(params.genes_file==""){
process GetGenesInfo{
   memory { strmem(params.other_mem_req) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 10
   output :
      file(out) into genes_file_ch
   publishDir "${params.output_dir}/data/",  mode:'copy'
   script :
     out="gencode.v19.genes"
     """
     wget -c ${params.genes_file_ftp} --no-check-certificate
     zcat `basename ${params.genes_file_ftp}` > file_genes
     change_genes_gencode.py file_genes
     """
}
}else{
genes_file_ch=Channel.fromPath(params.genes_file)
}



process MergeResult{
    label 'R'
    memory params.other_mem_req
    input :
      file(paintor) from res_paintor_ch
      file(infopaintor) from infores_paintor_ch
      file(cojo) from res_cojo
      file(caviarbf) from res_caviarbf
      tuple path(fmsss), path(fmssscred) from res_fmsss
      tuple path(fmcond), path(fmcondcred)  from res_fmcond
      file(datai) from data_i
      file(genes) from  genes_file_ch
      file(gwascat) from gwascat_ch
      file(pfileannot) from paintor_fileannotplot
   publishDir "${params.output_dir}/",  mode:'copy'
    output :
       set file("${out}.pdf"), file("${out}.all.out"), file("${out}.all.out")
    script :
      out=params.output
      infopaint=infopaintor.join(" ")
      //pfileannot= (baliseannotpaint==0) ? "":" --paintor_fileannot $pfileannot "
      pfileannot= " --paintor_fileannot $pfileannot "
      """
       cat $infopaint > infopaint
       echo "sss $fmsss" > infofinemap 
       echo "cond $fmcond" >> infofinemap 
       echo "cond $fmcondcred" > infofinemapcred
       echo "sss $fmssscred" >> infofinemapcred
       merge_finemapping_v2.r --out $out --listpaintor  infopaint  --cojo  $cojo --datai  $datai --caviarbf $caviarbf --list_genes $genes  --gwascat $gwascat --headbp_gc ${headgc_bp} --headchr_gc ${headgc_chr}  --listfinemap infofinemap  $pfileannot --listfinemap_cred infofinemapcred
      """
}

