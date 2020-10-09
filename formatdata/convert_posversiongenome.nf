
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      jean-tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2019
 *
 *
 * Description  : Nextflow pipeline to transform vcf file in plink and other format
 *
 *(C) University of the Witwatersrand, Johannesburg, 2016-2019 on behalf of the H3ABioNet Consortium
 *This is licensed under the MIT Licence. See the "LICENSE" file for details
 */


import java.nio.file.Paths;
import sun.nio.fs.UnixPath;
import java.security.MessageDigest;


// Checks if the file exists
allowed_params = ['file_toconvert','file_ref_gzip', "output_dir","output", "input_dir", "input_pat"]

/*file to convert if*/
params.file_toconvert=""
params.link_gwas_cat="https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
params.head_rs="SNPS"
params.head_chro="CHR_ID"
params.head_pos="CHR_POS"
params.output_dir='gwascat'
params.sep="TAB"
params.rs_info=""
params.poshead_rs_inforef=2
params.poshead_bp_inforef=3
params.file_ref_gzip=""
params.link_rs_info="ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz"
params.bin_crossmap="~/.local/bin/CrossMap.py"
params.data_crossmap=''
params.link_data_crossmap='http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'


if(params.file_toconvert==""){
process DlGwasCT{
   publishDir "${params.output_dir}/datai/", overwrite:true, mode:'copy'
   output :
       file("${fileconvert}err")
       file(fileconvert) into (file_convert_ch_ext, file_convert_ch)
   script :
      fileconverti="gwas_catalog_i.tsv"
      fileconvert="gwas_catalog.tsv"
      """
      wget -O $fileconverti ${params.link_gwas_cat}
      check_colfile.py $fileconverti $fileconvert
      """
}
}else{
file_convert_ch=Channel.fromPath(params.file_toconvert)
}
if(params.file_ref_gzip==""){
process DlInfoRs{
   output :
       file(file_rsinfo) into (file_rsinfo_ch)
   publishDir "${params.output_dir}/datai/", overwrite:true, mode:'copy'
   script :
      file_rsinfo="All_rs.vcf.gz"
      """
      wget -O $file_rsinfo ${params.link_rs_info}
      """
}

}else{
file_rsinfo_ch=Channel.fromPath(params.file_ref_gzip)
}

process ExtractInfo{
   input :
       file(fileconvert) from file_convert_ch_ext
   output :
       file("${headout}.rs") into rs_convert 
       file("${headout}.pos") into pos_convert
   script  :
      headout="search"
      """
      cp_extractpos.r --file $fileconvert  --out $headout --head_rs ${params.head_rs} --head_bp ${params.head_pos} --head_chr ${params.head_chro} --sep ${params.sep}
      """
}

process SearchPosWithRs{
   input :
     file(rsinfo) from file_rsinfo_ch
     file(rs_convert) from rs_convert
   publishDir "${params.output_dir}/tmpi/", overwrite:true, mode:'copy'
   output :
     file(outinfors) into outinfors_ch
   script :
     outinfors='info_extract.info'
     """
     zcat $rsinfo | cp_searchposwithrs.py  $rs_convert  ${outinfors}.tmp ${params.poshead_rs_inforef}
     awk '{if(NF>6){for(Cmt=7;Cmt<=NF;Cmt++)\$6=\$6";"\$Cmt};print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6}' ${outinfors}.tmp > $outinfors
     """
}

if(params.data_crossmap==""){
process DlDataCrossMap{
   publishDir "${params.output_dir}/datai/", overwrite:true, mode:'copy'
   output :
     file(fileout) into CrossMap_data_ch
   script :
    fileout=params.link_data_crossmap.split('/').last()
    """
    wget -c ${params.link_data_crossmap} 
    """
}
}else{
CrossMap_data_ch=Channel.fromPath(params.data_crossmap)
}

process CrossMapLaunch{
   input :
      file(CrossMapRef) from CrossMap_data_ch
      file(posI) from pos_convert
   publishDir "${params.output_dir}/tmpi/", overwrite:true, mode:'copy'
   output :
     file(poscrossmap) into result_crossmap
   script : 
    poscrossmap='convert_crossmap.cross'
    """
    ${params.bin_crossmap} bed $CrossMapRef $posI  $poscrossmap".tmp"
    awk '{if(NF>5){for(Cmt=6;Cmt<=NF;Cmt++)\$5=\$5";"\$Cmt};print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5}' $poscrossmap".tmp" > $poscrossmap
    """
}

process MergeRes{
   input :
    file(crossmap) from result_crossmap
    file(outinfors) from outinfors_ch
    file(filetoconvert) from file_convert_ch
   publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
   output :
     file("$headout*") 
   script : 
     headout=params.output
     """
     cp_mergepos.r --file $filetoconvert  --out $headout --head_rs ${params.head_rs} --head_bp ${params.head_pos} --head_chr ${params.head_chro} --sep ${params.sep} --file_rsres $outinfors --file_cross $crossmap
     """
}
