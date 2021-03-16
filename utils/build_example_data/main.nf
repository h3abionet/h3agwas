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
params.output="data"
params.output_dir="allgeno_data/"
params.pos_allgeno=""
params.nb_cpus=10
params.maf = 0.05
params.thin = 0.01
params.simu_hsq=0.3
params.simu_rep=10
params.gwas_cat=""
params.gwas_cat_ftp="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz"
params.list_pheno="Type 2 diabetes"
params.gcta_bin="gcta64"

params.sex_info_file="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
params.sex_info="gender,1:male,2:female,IID:sample"
params.simu_k=0.1
params.simu_cc_p=0.5
params.sex_error=0.05


params.clump_p1=0.0001 
params.clump_p2=0.01 
params.clump_r2=0.50 
params.clump_kb=250
params.dir_vcf="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"


if(params.list_chro_pheno==""){
params.list_chro_pheno=params.list_chro
}
listchro=getlistchro(params.list_chro)
listchro_pheno=getlistchro(params.list_chro_pheno)

listchro_ch=Channel.from(listchro)
listchro_ch2=Channel.from(listchro)
//Pattern100G="ALL.chr${chro}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"


listchro_ch=listchro_ch.combine(Channel.fromPath(params.pos_allgeno, checkIfExists:true))




process Dl1000G{
   cpus params.nb_cpus
   input :
      tuple val(chro), file(pos_geno) from listchro_ch 
   output :
       tuple val(chro), file("${file1000G}") into file1000G
   script :
      file1000G= (chro=='X') ? "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz" : "ALL.chr${chro}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      //ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
      file1000G= (chro=='Y') ? "ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz" : "$file1000G"
      """
      awk -v chro=$chro '{if(\$1==chro)print \$1"\\t"\$2-1"\\t"\$2"\\t"\$1":"\$2}' $pos_geno > pos.bed 
      bcftools view --threads ${params.nb_cpus} -R pos.bed ${params.dir_vcf}/$file1000G|bgzip -c > $file1000G
      """
}

process transfvcfInBed1000G{
   cpus params.nb_cpus
   input :
     tuple val(chro), file(vcf1000) from file1000G
   output :
     tuple val(chro), file("${out}.bim"), file("${out}.fam"), file("${out}.bed") into plk_chro
   script :
     out="chrtmp_"+chro
     """
     plink --vcf $vcf1000  --keep-allele-order  --make-bed -out $out --threads ${params.nb_cpus}
     """
}
//
process cleanPlinkFile{
    cpus params.nb_cpus
    input :
     tuple val(chro), file(bim), file(fam), file(bed) from plk_chro
    output : 
     tuple file("${out}.bim"), file("${out}.fam"), file("${out}.bed") into plk_chro_cl
    script :
     plk=bim.baseName
     out="chr_"+chro
     """
     cp "$bim" "${bim}.save"
     awk '{if(\$2=="."){\$2=\$1":"\$4};print \$0}' "${bim}.save" > "$bim"
     awk '{if(length(\$5)==1 && length(\$6)==1)print \$2}' $bim > ${bim}.wellpos.pos
     plink -bfile $plk  --keep-allele-order --extract  ${bim}.wellpos.pos --make-bed -out $out --threads ${params.nb_cpus}
     """
}
plk_chro_flt=plk_chro_cl.collect()

process mergePlinkFile{
   cpus params.nb_cpus
   input :
      val(allfile) from plk_chro_flt
   output :
     tuple file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into  allplkres_ch_befsex
   script : 
     println allfile
     allfile2=allfile.toList().collect {it.toString().replaceFirst(~/\.[^\.]+$/, '')}
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
     plink --bfile $firstbed  --keep-allele-order --threads ${params.nb_cpus} --merge-list listbedfile --make-bed --out $out
     fi

     """
}

process addSexFile{
 cpus params.nb_cpus
 input :
   tuple file(bed), file(bim), file(fam) from allplkres_ch_befsex
 publishDir "${params.output_dir}/geno_all/",  overwrite:true, mode:'copy'
 output :
    tuple file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into  allplkres_ch
    file("$out*")
 script :
     out=params.output
     basename=bed.baseName
     """
     wget -c ${params.sex_info_file}
     format_sex_forplk.r --file ${params.sex_info_file} --info ${params.sex_info} --out ${out}_sex --error ${params.sex_error}
     plink --bfile $basename --keep-allele-order --threads ${params.nb_cpus} --make-bed --out $out  --update-sex ${out}_sex
     """
}
process GwasCatDl{
    publishDir "${params.output_dir}/gwascat",  overwrite:true, mode:'copy'
    output :
       file("${out}.bed") into gwascat_bed
       file("${out}.pos") into gwascat_pos
       file("${out}_resume.csv") into gwascat_detail
       file("${out}*")
    script :
      out="gwascat_format"
      """
      wget -c ${params.gwas_cat_ftp} 
      format_gwascat.r --file `basename ${params.gwas_cat_ftp}` --pheno \"${params.list_pheno}\" --out $out  --chro ${listchro.join(',')}
      """
}

listchro_ch_gwascat=Channel.from(listchro_pheno)
listchro_ch_gwascat=listchro_ch_gwascat.combine(gwascat_pos)

process Dl1000G_GC{
   input :
      tuple val(chro), file(pos_geno) from listchro_ch_gwascat
   output :
       tuple val(chro), file("${file1000G}")  into file1000G_gwasc
   script :
      file1000G= (chro=='X') ? "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz" : "ALL.chr${chro}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      file1000G= (chro=='Y' ) ? "ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz" : "$file1000G"
      """
      tabix -fh ${params.dir_vcf}/$file1000G -R $pos_geno |bgzip -c > $file1000G
      """
}

process transfvcfInBed1000G_GC{
   cpus params.nb_cpus
   errorStrategy 'ignore'
   input :
     tuple val(chro), file(vcf1000)  from file1000G_gwasc
   output :
     tuple val(chro), file("${out}.bim"), file("${out}.fam"), file("${out}.bed") into plk_chro_gc
   script :
     out="chrtmp_"+chro
     """
     plink --vcf $vcf1000  --keep-allele-order  --make-bed -out $out --threads ${params.nb_cpus}
     """
}

process cleanPlinkFile_GC{
    cpus params.nb_cpus
    input :
     tuple val(chro), file(bim), file(fam), file(bed) from plk_chro_gc
    output : 
     tuple file("${out}.bim"), file("${out}.fam"), file("${out}.bed") into plk_chro_cl_gc
    script :
     plk=bim.baseName
     out="chr_"+chro
     """
     cp "$bim" "${bim}.save"
     awk '{if(\$2=="."){\$2=\$1":"\$4};print \$0}' "${bim}.save" > "$bim"
     awk '{if(length(\$5)==1 && length(\$6)==1)print \$2}' $bim > ${bim}.wellpos.pos 
     plink -bfile $plk  --keep-allele-order --extract  ${bim}.wellpos.pos --make-bed -out $out --threads ${params.nb_cpus}
     """
}


plk_chro_flt_gc=plk_chro_cl_gc.collect()

process mergePlinkFile_GC{
   cpus params.nb_cpus
   input :
      val(allfile) from plk_chro_flt_gc
   publishDir "${params.output_dir}/simul_pheno/datai/",  overwrite:true, mode:'copy'
   output :
     tuple file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into  allplkres_ch_gc
   script :
     println allfile
     allfile2=allfile.toList().collect {it.toString().replaceFirst(~/\.[^\.]+$/, '')}
     allfile2=allfile2.unique()
     firstbed=allfile2[0]
     allfile2.remove(0)
     nbfile=allfile2.size()
     out=params.output
     """
     if [ $nbfile == "0" ]
     then
     cp ${firstbed}.fam ${out}.fam
     cp ${firstbed}.bed ${out}.bed
     cp ${firstbed}.bim ${out}.bim
     else
     echo "${allfile2.join('\n')}" > listbedfile
     plink --bfile $firstbed  --keep-allele-order --threads ${params.nb_cpus} --merge-list listbedfile --make-bed --out $out
     fi
     """
}

process format_simulated{
   cpus params.nb_cpus
   input :
     tuple file(bed), file(bim), file(fam) from allplkres_ch_gc
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
     format_simulated.r --bfile $plk --gc_her $gwascat --out $outeffect --clump_p1 ${params.clump_p1} --clump_p2 ${params.clump_p2} --clump_r2 ${params.clump_r2} --clump_kb ${params.clump_kb}
     """
}


process simulation_quantitatif{
   cpus params.nb_cpus
   input :
     tuple file(bed), file(bim), file(fam), file(outeffect) from info_sim_qt
   publishDir "${params.output_dir}/simul_pheno/quant_pheno/",  overwrite:true, mode:'copy'
   output :
      file("$out")
   script :
     out=params.output+"_qt.pheno"
     plk=bed.baseName
     """
     ${params.gcta_bin} --bfile $plk --simu-causal-loci $outeffect  --simu-qt --simu-hsq ${params.simu_hsq} --out sim --simu-rep ${params.simu_rep}   --simu-k ${params.simu_k}
     format_file_sim.r --file sim".phen" --out $out   
     """
}

process simulation_qualitatif{
   cpus params.nb_cpus
   input :
     tuple file(bed), file(bim), file(fam), file(outeffect) from info_sim_ql
   publishDir "${params.output_dir}/simul_pheno/qual_pheno/",  overwrite:true, mode:'copy'
   output :
      file("$out")
   script :
     out=params.output+"_ql.pheno"
     plk=bed.baseName
     """
     ${params.gcta_bin} --bfile $plk --simu-causal-loci $outeffect  --simu-hsq ${params.simu_hsq} --out sim --simu-rep ${params.simu_rep}   --simu-k ${params.simu_k}  --simu-cc `estimated_cc.py $fam ${params.simu_k}`
     format_file_sim.r --file sim".phen" --out $out   
     """
}

