/* (C) H3ABionet
 GPL

 Scott Hazelhust, 2017-2018
*/


params.chipdescription = "/external/diskA/novartis/metadata/HumanOmni5-4v1_B_Physical-and-Genetic-Coordinates.txt"
inpat = "${params.input_dir}/${params.input_pat}"
params.output         = "chip"
params.strandreport   = false
params.manifest       = false




plink_src = Channel.create()


null_values = [0,"0",false,""]

if (! null_values.contains(params.samplesheet))
  samplesheet = Channel.fromPath(params.samplesheet)

def condChannel = { parm, descr ->
  filename = "/tmp/emptyZ0${descr}.txt";
  new File(filename).createNewFile()  
  if ((parm==0) || (parm=="0") || (parm == false) || (parm=="") || (parm=="false")) {
    return Channel.fromPath(filename)
  }  else {
    return  Channel.fromPath(parm);
  }
}


ref_ch      = condChannel(params.reference,"ref")
manifest_src_ch = condChannel(params.manifest,"man")
strand_src_ch   = condChannel(params.strandreport,"strn")

manifest_ch = Channel.create()
manifest1_ch = Channel.create()
manifest_src_ch.separate(manifest_ch,manifest1_ch)

strand_ch = Channel.create()
strand1_ch = Channel.create()
strand_src_ch.separate(strand_ch,strand1_ch) {a -> [a,a]}

ref1_ch = Channel.create()
ref2_ch = Channel.create()

ref_ch.separate(ref1_ch, ref2_ch) { a -> [a, a] }



def gChrom= { x ->
    println x
    m = x =~ /.*-(\w+).*/;
    return m[0][1]
}



  output      = params.output
  // array may be the manifest or pref a file with both genetic 
  // and physical coordinates
  array       = Channel.fromPath(params.chipdescription)
  report       = Channel.fromPath(inpat)
  output_align = params.output_align




  

  process illumina2lgen {
    maxForks params.max_forks
    input:
       set file(report), file(array) from report.combine(array)
    output:
       set file("${output}.ped"), file("${output}.map") into ped_ch
    script:
        idpat = params.idpat
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
     set file("raw.bed"), file("raw.fam"), file("raw.log") into plink_src
     set file("rawraw.bim") into fill_in_bim_ch
     file ('raw.fam') into fix_fam_ch
   script:
    """
    ls *.bed | sort > beds
    ls *.bim | sort >  bims
    ls *.fam | sort >  fams
    paste beds bims fams >  mergelist
    plink --merge-list mergelist --make-bed --out raw
    mv raw.bim rawraw.bim
    """
 } 

  process fillInBim {  //  Deals with monomorphic or non-called SNPs
    input:
     file(inbim) from fill_in_bim_ch
     file(strand) from strand1_ch
     file(manifest) from manifest1_ch
    output:
     file("raw.bim") into filled_bim_ch
    script:
       "fill_in_bim.py $strand $manifest $inbim  raw.bim"
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



 process checkReferenceFormat {

   input:
     file(ref) from ref1_ch
   output:
     stdout       into ref_parm_ch
   script:
      refBase = ref.baseName
      newref = "${refBase}.new"
      "checkRef.py  $ref $newref"
   }


 process alignStrand {
   input:
    set file(bed), file(fam), file(logfile) from plink_src
    file(bim) from filled_bim_ch
    file(ref) from ref2_ch
    val(ref_parm) from ref_parm_ch
    file(flips) from flip_ch
    publishDir params.output_dir, pattern: "*.{bed,bim,log,badsnps}", \
        overwrite:true, mode:'copy'
   output:
    set file("*.{bed,bim,log,badsnps}") into aligned_ch
   script:
    base = bed.baseName
    refBase = ref.baseName
    ref_parm = ref_parm.trim()
    opt = "--a2-allele $ref $ref_parm"
    if (refBase=="empty") opt="--keep-allele-order"
    """
    plink --bfile $base $opt --flip $flips --make-bed --out $output
    grep Impossible ${output}.log | tr -d . | sed 's/.*variant//'  > ${output}.badsnps
    """
 }




if (null_values.contains(params.samplesheet)) {
  process fixFam{
    input:
    file(fam) from fix_fam_ch
   publishDir params.output_dir, pattern: "${output}.fam", \
             overwrite:true, mode:'copy'
   output:
     set file("${output}.fam") into fixedfam_ch
   script:
     "cp $fam ${output}.fam"
  }
} else {
 process fixFam {
   input:
     file(samplesheet)
     file(fam) from fix_fam_ch
   publishDir params.output_dir, pattern: "${output}.fam", \
             overwrite:true, mode:'copy'
   output:
     set file("${output}.fam") into fixedfam_ch
   script:
     idpat = params.idpat
     batch_col = params.batch_col
     template "sheet2fam.py"
 } 
}

fixedfam_ch.combine(aligned_ch).subscribe { 
  files =  it.collect { fn -> fn.getName().replace(".*/","") }
  println "The output can be found in ${params.output_dir}: files are ${files}"
}


