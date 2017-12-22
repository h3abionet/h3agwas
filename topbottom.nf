/* (C) H3ABionet
 GPL

 Don Armstrong
 Scott Hazelhrust
*/

params.vcf_merge_table="/data/aux/RsMergeArch.bcp.gz"
params.reference      ="/data/aux/Homo_sapiens.GRCh37.75.dna.toplevel.fa"
params.dbsnp_all      ="/data/aux/All_20170403.vcf.gz"
params.manifest       ="/data/aux/HumanOmni5-4-v1-0-D.csv"
params.chipdescription = "/external/diskA/novartis/metadata/HumanOmni5-4v1_B_Physical-and-Genetic-Coordinates.txt"
inpat = "${params.input_dir}/${params.input_pat}"
params.output         = "chip"


params.vcf = false

plink_src = Channel.create()

checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}



def gChrom= { x ->
    println x
    m = x =~ /.*-(\w+).*/;
    return m[0][1]
}



// Take a list of channels and phase them
def combineElts = { u, v ->
  ut=u.getClass()
  vt=v.getClass()
  ft = file("tmp").getClass()
  lt = [].getClass()
  st = "tmp".getClass()
  inty = 1.getClass()
  // if the two arguments are singletons return the list
  if  ( (ut in [st,inty,ft]) && vt in [st,inty,ft]) {
    return [u,v]
  }
  else
    if (ut != lt)   u = [u]
    if (vt == lt) 
       if (v[0].getClass() == ft) 
         return u+v;
       else 
         // The zeroth element of v is a key and we don't want to repeat
         return u+v[1..-1]
    else 
    return u+[v];
 
}

def phaseAllLists = { channels ->
     start = channels[0]
     println("HELLO")
     for (c in channels[1..-1]){
         start = start.phase(c) .map(combineElts)
     }
     return start;
}




  output      = params.output

  fam         = Channel.fromPath(params.fam)

  rsmerge_ch  = Channel.create()


  ref_ch      = Channel.fromPath(params.reference)
  remap       = Channel.fromPath("scripts/remap_map_and_ped.pl")
  dbsnp_all   = Channel.fromPath(params.dbsnp_all)
  array        = Channel.fromPath(params.chipdescription)
  array2        = Channel.fromPath(params.chipdescription)
  manifest_csv   = Channel.fromPath(params.manifest)
  report        = Channel.fromPath(inpat)
  vcf_rs_grep  = file("templates/vcf_rs_grep")
  remap        = file("templates/remap_map_and_ped.pl")
  fix_strandedness = Channel.fromPath("scripts/fix_strandedness.pl")
  chroms = 1..26
  (1..26).each {   rsmerge_ch.bind(file(params.vcf_merge_table)) }



  process illumina2lgen {
    input:
       file(report)
       file(array)
       file(fam)
    output:
       file("$output*") into lgenfs
    script:
        report = inpat
        template "topbot2plink.py"
  }




  process lgen2ped {
    
   input:
     set val(base), file (fs) from lgenfs.flatten().map { fn -> tuple (gChrom(fn) , fn) }.groupTuple ()
   output:
     set val(base), file("${base}.ped"), file("${base}.map") into plink_chroms
     set val(base), file("${base}.map") into  map_ch1, map_ch2
   script:
      """
      echo $base
      echo $fs

      plink --lgen chip-${base}.lgen --map chip-${base}.map --fam chip-${base}.fam --recode --out $base
      """
  }


  process getChromSnps {
    // extract out SNPs for each chromosome
    time '10h' 
    input:
       file(dbsnp_all) 
    output:
       file("*-*.vcf.gz") into dbsnp_chrom
    script:
      inputf=dbsnp_all
      template "vcf_split_chrom.py"
      //"""
      //for x in `seq 1 22`; do
      //    touch All_20170403-\${x}.vcf.gz
      //done
      //touch All_20170403-X.vcf.gz
      //touch All_20170403-Y.vcf.gz
      //touch All_20170403-MT.vcf.gz
      //"""
  }

  process grepAnnotation { 
    input:
      set val(base), file(mapf) from map_ch1
    output:
      set val(base), file("${base}_grep") into base_grep
    script:
         "awk '{print \$\$2}' $mapf > ${base}_grep"
  }


  process annotation {
    input:
      set val(chrom), file(dbsnp), file(awkout) from \
            phaseAllLists([dbsnp_chrom.flatten().map { b -> [gChrom(b), b]},base_grep])
      file (vcf_merge_table) from  rsmerge_ch
    output:
      set val(chrom), file("annotation.vcf.gz"), file(vcf_merge_table) into annotation_ch
    script:
      chrom_vcf=dbsnp_chrom
      "${vcf_rs_grep}  $vcf_merge_table $dbsnp $awkout | gzip -c > annotation.vcf.gz"
    }


  process fix_plink {
    input:
      set val(base), file(ped), file(map), file(annotation), file(bcp) from phaseAllLists([plink_chroms,annotation_ch])
    output:
      set val(base), file("${base}_fixed.ped"), file("${base}_fixed.map"), file(annotation), file(bcp) into fixed_ch
    script:
    base = ped.baseName    
    """
       $remap --vcf $annotation --ped $ped --map ${map} \
              --ped_out ${base}_fixed.ped --map_out ${base}_fixed.map --merge $bcp
    """
  }


  process flip {
    input:
    set val(base), file(ped), file(map), file(annotation), file(bcp), \
        file(fix_strandedness), file(array2), file(ref_genome) from \
            fixed_ch.combine(fix_strandedness).combine(array2).combine(ref_ch)
    output:  
       set file("${flipped}.ped"),file("${flipped}.map"),file("${flipped}_flipping.log") into create_bed_ch
    publishDir params.output_dir, pattern: "*log", \
             overwrite:true, mode:'copy'
    script:
       flipped = "${base}_fl"
       """
	   ./${fix_strandedness} --ped $ped --map $map                \
	   --ped-out ${flipped}.ped --map-out ${flipped}.map		\
           --illumina-data ${array2} \
           --merge $bcp  \
           --vcf $annotation \
           --ref $ref_genome                                       \
           > ${flipped}_flipping.log          \
      """
    }	  

 process create_bed {
   input:
     file(data) from create_bed_ch.toList()
   output:
     set file("out.{bed,bim,fam}") into plink_src
   script:
    """
    ls *.ped | sort > peds
    ls *.map | sort >  maps
    paste peds maps >  mergelist
    plink --merge-list mergelist --make-bed --out $output
    """
 } 


