/* (C) H3ABionet
 GPL

 Scott Hazelhust, 2017-2018
*/


params.chipdescription = "/external/diskA/novartis/metadata/HumanOmni5-4v1_B_Physical-and-Genetic-Coordinates.txt"
inpat = "${params.input_dir}/${params.input_pat}"
params.output         = "chip"
params.strandreport   = false
params.manifest       = false
params.idpat         = "(.*)"



abbrev = params.abbrev

plink_src = Channel.create()

def condChannel = { parm, descr ->
  if ((parm==0) || (parm=="0") || (parm == false) || (parm=="false")) {
    println(parm+" "+descr)
    File empty = new File("empty."+descr)
    empty.createNewFile()
    return Channel.fromPath(empty)
  }  else {
    return  Channel.fromPath(parm);
  }
}


ref_ch      = condChannel(params.reference,"ref")
manifest_ch = condChannel(params.manifest,"man")
strand_ch   = condChannel(params.strandreport,"strn")





def gChrom= { x ->
    println x
    m = x =~ /.*-(\w+).*/;
    return m[0][1]
}



  output      = params.output
  samplesheet = Channel.fromPath(params.samplesheet)
  // array may be the manifest or pref a file with both genetic 
  // and physical coordinates
  array       = Channel.fromPath(params.chipdescription)
  report       = Channel.fromPath(inpat)
  output_align = params.output_align




  

  process illumina2lgen {
    maxForks 95
    input:
       set file(report), file(array) from report.combine(array)
    output:
       set file("${output}.ped"), file("${output}.map") into ped_ch
    script:
        samplesize = params.samplesize
	idpat      = params.idpat
        output = report.baseName
        template "topbot2plink.py"
  }




  process bedify {
    
   input:
     set file(ped), file(map) from ped_ch
   output:
     file("${base}.bed") into bed_ch
     file ("${base}.bim") into bim_ch
     file("${base}.fam") into fam_ch
   script:
      base = ped.baseName
      """
      echo $base

      plink --file $base --make-bed --out $base
      """
  }



 process combineBed {
   input:
     file(bed) from bed_ch.toList()
     file(bim) from bim_ch.toList()
     file(fam) from fam_ch.toList()
   output:
     set file("raw.bed"), file("raw.bim"), file("raw.fam"), file("raw.log") into plink_src
     file ('raw.fam') into fix_fam_ch
   script:
    """
    ls *.bed | sort > beds
    ls *.bim | sort >  bims
    ls *.fam | sort >  fams
    paste beds bims fams >  mergelist
    plink --merge-list mergelist --make-bed --out raw
    """
 } 


 process getFlips {
   input:
     file(strandreport) from strand_ch
     file(manifest) from manifest_ch
   output:
     file(flips) into flip_ch
   script:
    flips = "flips.lst"
    template "strandmismatch.py"
 }


 process alignStrand {
   input:
    set file(bed), file(bim), file(fam), file(logfile) from plink_src
    file(ref) from ref_ch
    file(flips) from flip_ch
    publishDir params.output_dir, pattern: "*.{bed,bim,log,badsnps}", \
        overwrite:true, mode:'copy'
   output:
    set file("*.{bed,bim,log,badsnps}") into aligned_ch
   script:
    base = bed.baseName
    refBase = ref.baseName
    opt = "--a2-allele $ref"
    if (refBase=="empty") opt="--keep-allele-order"
    """
    plink --bfile $base $opt --flip $flips --make-bed --out $output
    grep Impossible ${output}.log | tr -d . | sed 's/.*variant//'  > ${output}.badsnps
    """
 }


 process fixFam {
   input:
     file(samplesheet)
     file(fam) from fix_fam_ch
   publishDir params.output_dir, pattern: "${output}.fam", \
             overwrite:true, mode:'copy'
   output:
     set file("${output}.fam") into fixedfam_ch
   script:
    batch_col = params.batch_col
    template "sheet2fam.py"
 } 


fixedfam_ch.combine(aligned_ch).subscribe { 
  files =  it.collect { fn -> fn.getName().replace(".*/","") }
  println "The output can be found in ${params.output_dir}: files are ${files}"
}


