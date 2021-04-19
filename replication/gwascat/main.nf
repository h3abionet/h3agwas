#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths



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


def helps = [ 'help' : 'help' ]
allowed_params = ["cut_maf", "output_dir", "pb_around_rs", "mem_req", "work_dir","mem_req","big_time", "output","nb_cpu" , "input_dir","input_pat", "file_gwas", "gwas_cat", "site_wind", "r2_clump","min_pval_clump", "size_win_kb"]
allowed_params_blocks = ["haploblocks", "plkref_haploblocks", "plk_othopt_haploblocks"]
allowed_params_other=["max_forks", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key","region", "AMI","maxInstances","instance-type", "instanceType", "bootStorageSize", "boot-storage-size", "max-instances", "sharedStorageMount", "shared-storage-mount", "scripts"]
allowed_params_headinfo=["head_chro_gwascat", "head_bp_gwascat", "head_pval_gwascat"]
allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A0"]
allowed_params+=allowed_params_head
allowed_params+=allowed_params_other
allowed_params+=allowed_params_blocks
allowed_params+=allowed_params_headinfo
params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}

checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}



def params_help = new LinkedHashMap(helps)
params.queue      = 'batch'
params.work_dir   = "$PWD"
params.input_listfiles = "${params.work_dir}/list_files.input"
params.output_dir = "${params.work_dir}/output"
params.output = "replication"
params.cut_maf = 0.01


params.mem_req="8G"
params.big_time="1000H"

params.head_pval = "P_BOLT_LMM"
params.head_freq = ""
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.data = ""
params.head_beta=""
params.head_se=""
params.head_A1="ALLELE1"
params.head_A0="ALLELE0"

#params.clump_p1=0.01
#params.clump_p2=1
params.clump_r2=0.1


params.max_pval_rep=10**-6
params.size_win_kb=25
params.nb_cpu = 3

params.gwas_cat=""
// haploblocks information
params.haploblocks=""
params.plkref_haploblocks=""
params.plk_othopt_haploblocks=""

//params plink
params.plink_bin='plink'
params.genes_file=""
params.max_plink_cores = 4

max_plink_cores=params.max_plink_cores 
plink_mem_req = params.plink_mem_req


//params.gene_file_ftp="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"

//params gwas cat 
params.head_bp_gwascat="bp"
params.head_chro_gwascat="chro"
params.head_pval_gwascat="pvalue"
params.head_af_gwascat="risk.allele.af"
params.head_rs_gwascat="rs"
params.gwas_cat_ftp="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz"
params.list_pheno="Type 2 diabetes"
params.list_chro="1-22"
if(params.input_dir=="" || params.input_pat==""){
println "params input_dir directory of your bedfile or/and input_pat pattern of your bedfile not define"
System.exit(-1)

}

if(params.head_beta==""){
println "beta header of your gwas file must be provide : --head_beta"
System.exit(-1)
}
if(params.head_se==""){
println "se header of your gwas file must be provide : --head_se"
System.exit(-1)
}

if(params.head_rs==""){
println "rs header of your gwas file must be provide : --head_rs"
System.exit(-1)
}

if(params.head_bp==""){
println "pos header of your gwas file must be provide : --head_bp"
System.exit(-1)
}



//params.info_gwascat="DISEASE.TRAIT,REPORTED.GENE.S.,MAPPED_GENE,INITIAL.SAMPLE.SIZE"
params.threshold_pval_gwascat=1
params.threshpval=0.05

 
//BedFileI=/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/locuszoom/data/1000G/genotypes/2014-10-14/EUR/chr$Chro
//plink -bfile $BedFileI --blocks no-pheno-req --out $DirOut/chr$Chro"_block"  &

params.r2_clump=0.1
params.min_pval_clump=0.1

if(params.gwas_cat==""){
listchro=getlistchro(params.list_chro)
process get_gwascat_hg19{
   label 'R'
   publishDir "${params.output_dir}/gwascat",  overwrite:true, mode:'copy'
   output :
       file("${out}.bed") into gwascat_bed
       file("${out}_range.bed") into gwascat_rangebed
       file("${out}.pos") into gwascat_pos
       file("${out}_resume.csv") into gwascat_detail_statpos
       file("${out}_all.csv") into (gwascat_all_statpos,gwascat_all_statld, gwascat_all_statclump,gwascat_all_statwind)
       file("${out}*")
   script :
     chroparam= (params.list_chro=='') ? "" : " --chro ${listchro.join(',')}"
     phenoparam= (params.list_pheno=='') ? "" : " --pheno \"${params.list_pheno}\" "
     out=params.output
   """
   wget -c ${params.gwas_cat_ftp}
   format_gwascat.r --file `basename ${params.gwas_cat_ftp}` $chroparam $phenoparam --out $out --wind ${params.size_win_kb} 
   """
}
gwascathead_chr="chrom"
gwascathead_bp="chromEnd"
infogwascat="pubMedID;author;trait;initSample"

}

filegwas_chrextr=Channel.fromPath(params.file_gwas,checkIfExists:true)


process extractgwas_fromgwascat{
   input : 
     file(pos) from gwascat_pos
     file(gwas) from filegwas_chrextr
   publishDir "${params.output_dir}/gwas_sub",  overwrite:true, mode:'copy'
   output :
     file("${params.output}_range.assoc") into (clump_file_ch,ld_file_ch)
     file("${params.output}_pos.init") into (pos_file_ch, range_file_ch)
    file("${params.output}_range.init") into (range_file_ch_clump,wind_file_ch)
     file("${params.output}*")
   script :
    """
    extract_posgwas.py --bed $pos --gwas $gwas --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --wind ${params.size_win_kb}  --pval_gwas ${params.head_pval} --rs_gwas ${params.head_rs}  --out ${params.output}
    """
}

bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()
raw_src_ch= Channel.create()

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set {raw_src_ch}

plk_ch_clump  = Channel.create()
plk_ch_ld  = Channel.create()
plk_ch_clumpstat  = Channel.create()

    /*JT : append boltlmm_assoc_ch and a]*/
raw_src_ch.separate(plk_ch_clump,plk_ch_ld, plk_ch_clumpstat) { a -> [a,a,a] }



process clump_aroundgwascat{
    cpus max_plink_cores
    memory plink_mem_req
   input :
      file(assocclump) from clump_file_ch
      set file(bed), file(bim), file(fam) from plk_ch_clump
   publishDir "${params.output_dir}/res/clump/tmp",  overwrite:true, mode:'copy'
   output :
      file("${out}.clumped") into clump_res_ch
      file("$out*")   
   script :
      bfile=bed.baseName
      out=params.output
      """ 
      plink -bfile $bfile  --clump $assocclump -clump-p1 $params.min_pval_clump --clump-p2 1 --clump-kb ${params.size_win_kb} --clump-r2 $params.clump_r2 -out $out --threads $max_plink_cores
      """
}

process computedstat_pos{
   label 'R'
   input :
        file(assocpos) from pos_file_ch
        file(gwascat)  from gwascat_all_statpos
   publishDir "${params.output_dir}/res/pos",  overwrite:true, mode:'copy'
   output :
      set file("${out}.csv"), file("${out}_cmpfrequencies.svg"), file("${out}_cmpz.svg") 
      file("$out*") 
   script :
    out=params.output+"_pos"
    af= (params.head_freq=='') ? "" : " --af_gwas ${params.head_freq} "
    """
    computestat_pos.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  $af --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out
    """
}


process computedstat_win{
   label 'R'
   input :
        file(assocpos) from wind_file_ch
        file(gwascat)  from gwascat_all_statwind
   publishDir "${params.output_dir}/res/wind",  overwrite:true, mode:'copy'
   output :
      set file("${out}.csv"), file("${out}_cmpfrequencies.svg"), file("${out}_cmpz.svg")
      file("$out*")
   script :
    out=params.output+"_wind"
    af= (params.head_freq=='') ? "" : " --af_gwas ${params.head_freq} "
    """
    computestat_wind.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  $af --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --min_pval ${params.threshold_pval_gwascat} --info_gwascat  \"$infogwascat\"
    """
}



process computed_ld{
    cpus max_plink_cores
    memory plink_mem_req
  input :
      file(filegwascat) from gwascat_rangebed
      set file(bed), file(bim), file(fam) from plk_ch_ld
    publishDir "${params.output_dir}/res/ld/tmp",  overwrite:true, mode:'copy'
    output :
       file("${out}.ld") into ld_res_ch
       file("$out*") 
    script :
    out=params.output+"_ld"
    plkf=bed.baseName
    """
    plink -bfile $plkf  --r2 --extract range  ${filegwascat}  --ld-window-kb $params.size_win_kb        --ld-window-r2 $params.clump_r2 -out $out --threads $max_plink_cores
    """

}


process computed_ld_stat{
    label 'R'
    input :
      file(fileld) from  ld_res_ch
      file(gwascat) from gwascat_all_statld
      file(assocpos) from range_file_ch
    publishDir "${params.output_dir}/res/ld/",  overwrite:true, mode:'copy'
    output :
       file("$out*")
    script :      
      out=params.output+"_ld"
      af= (params.head_freq=='') ? "" : " --af_gwas ${params.head_freq} "
      """
      computestat_ld.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  $af --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --ld_file $fileld --min_pvalue ${params.min_pval_clump} --min_r2  ${params.clump_r2} --info_gwascat \"$infogwascat\"
      """
}

process computed_clump_stat{
    label 'R'
    input :
      file(fileclum) from  clump_res_ch
      file(gwascat) from gwascat_all_statclump
      file(assocpos) from range_file_ch_clump
      set file(bed), file(bim), file(fam) from plk_ch_clumpstat

    publishDir "${params.output_dir}/res/clump/",  overwrite:true, mode:'copy'
    output :
       file("$out*")
    script :
      out=params.output+"_ld"
      af= (params.head_freq=='') ? "" : " --af_gwas ${params.head_freq} "
      """
      computestat_clump.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  $af --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --clump_file $fileclum --min_pvalue ${params.min_pval_clump} --min_r2  ${params.clump_r2} --info_gwascat \"$infogwascat\" --bim $bim
      """

}
