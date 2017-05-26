/* (C) H3ABionet
 GPL

 Don Armstrong
 Scott Hazelhrust
*/

params.vcf_merge_table="/data/aux/RsMergeArch.bcp.gz"
params.reference      ="/data/aux/Homo_sapiens.GRCh37.75.dna.toplevel.fa"
params.dbsnp_all_vcf  ="/data/aux/All_20170403.vcf.gz"
params.manifest       ="/data/aux/HumanOmni5-4-v1-0-D.csv"

inpat = "${params.input_dir}/${params.input_pat}"



params.vcf = false

plink_src = Channel.create()

checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}


if (params.topbot) {

  topbot_chan = Channel.fromPath("${inpat}.lgen")
  lmap = Channel.fromPath("${inpat}.map")
  fam         = Channel.fromPath(params.fam)
  rsmerge_ch  = Channel.fromPath(params.vcf_merge_table)
  ref_ch      = Channel.fromPath(params.reference)
  vcf_rs_grep = Channel.fromPath("scripts/vcf_rs_grep")
  remap       = Channel.fromPath("scripts/remap_map_and_ped.pl")
  all_vcf     = Channel.fromPath(params.dbsnp_all_vcf)
  fix_strandedness = Channel.fromPath("scripts/fix_strandedness.pl")
  array_csv   = Channel.fromPath(params.manifest)



  process lgen2ped {
    
    input:
      file(topbot_chan)
      file(lmap)
      file(fam)
    output:
      set file("${base}.ped"), file("${base}.map") into input_ch
      file("${base}.map") into mapf
   script:
      base = topbot_chan.baseName
      """
      plink --lgen $topbot_chan --map $lmap --fam $fam --recode --out $base
      """
  }


  process annotation {
    input:
       file(mapf)
       file(vcf_rs_grep)
       file(vcf_merge_table) from rsmerge_ch
       file(all_vcf) 
    output:
    set file("annotation_vcf.gz"), file("${vcf_merge_table}") into annotation_ch
    script:
      """
       awk '{print \$\$2}' $mapf | \
       ./vcf_rs_grep --merge  $vcf_merge_table  $all_vcf |gzip -c > annotation_vcf.gz
      """
    }


  process fix_plink {
    input:
      set file(ped), file(map) from  input_ch
      file(remap)
      set file(annotation), file(bcp) from annotation_ch
    output:
      set val(base), file("${base}_fixed.ped"), file("${base}_fixed.map"), file(annotation), file(bcp) into fixed_ch
    script:
    base = ped.baseName    
    """
       ./$remap --vcf $annotation --ped $ped --map ${map} \
             --ped_out ${base}_fixed.ped --map_out ${base}_fixed.map --merge $bcp
    """
  }

  
  process flip {
    input:
       set val(base), file(ped), file(map), file(annotation), file(bcp) from fixed_ch
      file(fix_strandedness)
      file(array_csv)
      file(ref_genome) from ref_ch
    output:  
      set val(base), file("${flipped}.ped"),file("${flipped}.map"),file("${flipped}_flipping.log") into create_bed_ch
    publishDir params.output_dir, pattern: "*log", \
             overwrite:true, mode:'copy'
    script:
       flipped = "${base}_fl"
       """
	   ./${fix_strandedness} --ped $ped --map $map                \
	   --ped-out ${flipped}.ped --map-out ${flipped}.map		\
           --illumina-data ${array_csv} \
           --merge $bcp  \
           --vcf $annotation \
           --ref $ref_genome                                       \
           > ${flipped}_flipping.log          \
      """
    }	  

 process create_bed {
   input:
     set val(base), file(ped), file(map) from create_bed_ch
   output:
     set file("$base.{bed,bim,fam}") into plink_src
   script:
    """
    plink --bed $ped --fam $map --make-bed --out $base
    """
 } 

} else {

Channel
.fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName }  \
   .ifEmpty { error "No matching plink files" }        \
   .map { a -> [checker(a[1]), checker(a[2]), checker(a[3])] }
   .into (plink_src)  							\



} 





plink_src.subscribe { println "$it" }