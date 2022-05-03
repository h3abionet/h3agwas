#!/usr/bin/env nextflow
import java.nio.file.Paths
/*
 * Authors       :
 *
 *
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2019
 *
 *
 * Description  : Nextflow pipeline for extract datatuple of 1000 genomes
 *
 */

//---- General definitions --------------------------------------------------//
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


params.list_chro="1-22,X"
params.list_chro_pheno=""

params.simu_hsq=0.3
params.simu_rep=10
params.gcta_bin="gcta64"


params.gwas_cat=""
params.gwas_cat_ftp="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz"
params.list_pheno="Type 2 diabetes"

params.simu_k=0.1
params.simu_cc_p=0.5
params.sex_error=0.05


params.input_pat=""
params.input_dir=""

params.clump_p1=0.0001
params.clump_p2=0.01
params.clump_r2=0.50
params.clump_kb=250
params.nb_snp=-1

params.nb_cpus=4

if(params.list_chro_pheno==""){
list_chro_pheno=params.list_chro
}else{
list_chro_pheno=params.list_chro_pheno
}

listchro=getlistchro(params.list_chro)
listchro_pheno=getlistchro(list_chro_pheno)

listchro_ch=Channel.from(listchro)
listchro_ch2=Channel.from(listchro)


// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}



raw_src_ch= Channel.create()

bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }


gwas_plk_ch = Channel.create()
gcta_plk_ch = Channel.create()
raw_src_ch.separate(gwas_plk_ch, gcta_plk_ch) { a -> [a,a] }



process GwasCatDl{
    label 'R'
    publishDir "${params.output_dir}/gwascat",  overwrite:true, mode:'copy'
    output :
       file("${out}.bed") into (gwascat_bed_ch1,gwascat_bed_ch2)
       file("${out}.pos") into gwascat_pos
       file("${out}_resume.csv") into gwascat_detail
       file("${out}*")
    script :
      out="gwascat_format"
      """
      wget -c ${params.gwas_cat_ftp} --no-check-certificate
      format_gwascat.r --file `basename ${params.gwas_cat_ftp}` --pheno \"${params.list_pheno}\" --out $out  --chro ${listchro.join(',')}
      """
}


process extract_gc_dl{
   cpus params.nb_cpus
   input :
     tuple file(bed), file(bim), file(fam) from gwas_plk_ch
     file(bed_pos) from gwascat_bed_ch1
   output : 
     set file("${out}.bed"), file("${out}.bim"),  file("${out}.fam") into gwas_plk_ch_gc
   script :
      plk=bed.baseName
      out=plk+"_gc"
      """
      plink -bfile $plk --extract range $bed_pos --keep-allele-order  --make-bed -out $out
      """
}

process format_simulated{
   label 'R'
   cpus params.nb_cpus
   input :
     tuple file(bed), file(bim), file(fam) from gwas_plk_ch_gc
     file(gwascat) from gwascat_detail
   publishDir "${params.output_dir}/simul_pheno/pheno_format/",  overwrite:true, mode:'copy'
   output :
      tuple file(bed), file(bim), file(fam), file(outeffect) into (info_sim_qt, info_sim_ql)
      file(fam) into fam_countnb
      file("$params.output*")
   script :
     plk=bed.baseName
     outeffect=params.output+".effect.rs"
     """
     format_simulated.r --bfile $plk --gc_her $gwascat --out $outeffect --clump_p1 ${params.clump_p1} --clump_p2 ${params.clump_p2} --clump_r2 ${params.clump_r2} --clump_kb ${params.clump_kb} --nb_snp ${params.nb_snp}
     """
}

process simulation_quantitatif{
   label 'gcta'
   cpus params.nb_cpus
   input :
     tuple file(bed), file(bim), file(fam), file(outeffect) from info_sim_qt
   output :
      file("sim.phen") into sim_ql 
   script :
     out=params.output+"_qt.pheno"
     plk=bed.baseName
     """
     ${params.gcta_bin} --bfile $plk --simu-causal-loci $outeffect  --simu-qt --simu-hsq ${params.simu_hsq} --out sim --simu-rep ${params.simu_rep}   --simu-k ${params.simu_k}
     """
}
process format_sim_quantitatif{
    label 'R'
    input :
       file(file) from sim_ql
   publishDir "${params.output_dir}/simul_pheno/quant_pheno/",  overwrite:true, mode:'copy'
    output :
      file("$out")
    script :
      out=params.output+"_qt.pheno"
      """
      format_file_sim.r --file $file --out $out
      """
}

process simulation_qualitatif{
   label 'gcta'
   cpus params.nb_cpus
   input :
     tuple file(bed), file(bim), file(fam), file(outeffect) from info_sim_ql
   output :
      file("sim.phen") into sim_qt
   script :
     out=params.output+"_ql.pheno"
     plk=bed.baseName
     """
     ${params.gcta_bin} --bfile $plk --simu-causal-loci $outeffect  --simu-hsq ${params.simu_hsq} --out sim --simu-rep ${params.simu_rep}   --simu-k ${params.simu_k}  --simu-cc `estimated_cc.py $fam ${params.simu_k}`
     """
}


process format_sim_qualitatif{
    label 'R'
    input :
       file(file) from sim_qt
   publishDir "${params.output_dir}/simul_pheno/qual_pheno/",  overwrite:true, mode:'copy'
    output :
      file("$out")
    script :
      out=params.output+"_ql.pheno"
      """
      format_file_sim.r --file $file --out $out
      """
}

