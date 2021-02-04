/* 

(C) University of the Witwatersrand, Johannesburg, 2016-2018 on behalf of the H3ABioNet Consortium

 This is licensed under the MIT Licence. See the "LICENSE" file for details

 Scott Hazelhust, 2017-2018
*/


null_values = [0,"0",false,""]

params.chipdescription = "/external/diskA/novartis/metadata/HumanOmni5-4v1_B_Physical-and-Genetic-Coordinates.txt"
inpat = "${params.input_dir}/${params.input_pat}"
params.output         = "chip"
params.strandreport   = false
params.manifest       = false
params.queue          = 'batch'
params.output_align   = 'ref'
params.mask_type      = 0
params.sheet_columns  = 0

params.idpat=0
params.indiv_memory_req="2GB"
params.combined_memory_req="8GB"
params.time_req="12h"

params.newpat = 0
params.mask = 0
params.replicates = 0
mask = params.mask
idpat = params.idpat

mask_type = params.mask_type

if (!null_values.contains(params.sheet_columns))
  sheet_cols_ch = Channel.fromPath(params.sheet_columns)
else 
  sheet_cols_ch =  Channel.fromPath("0")

plink_src = Channel.create()




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
manifest_src_ch.separate(manifest_ch,manifest1_ch) { a-> [a,a] }

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
  array       =  Channel.fromPath(params.chipdescription)
  report       = Channel.fromPath(inpat).ifEmpty { error "No files match the pattern "+inpat }
  output_align = params.output_align




  

  process illumina2lgen {
    maxForks params.max_forks
    memory params.indiv_memory_req
    time   params.time_req
    input:
       set file(report), file(array) from report.combine(array)
    output:
       set file("${output}.ped"), file("${output}.map") into ped_ch
    script:
        samplesize = params.samplesize
	idpat      = params.idpat
        output = report.baseName
        """
        hostname
        topbottom.py $array $report $samplesize '$idpat' $output
        """
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

      plink --file $base --no-fid --make-bed --out $base
      """
  }



// Combine  all plinks and possibly remove errprs
 process combineBed {
   input:
     file(bed) from bed_ch.toList()
     file(bim) from bim_ch.toList()
     file(fam) from fam_ch.toList()
   memory params.combined_memory_req
   time   params.time_req
   output:
     set file("raw.bed"), file("raw.fam"), file("raw.log") into plink_src
     file("rawraw.bim") into fill_in_bim_ch
   script:
    """
    ls *.bed | sort > beds
    ls *.bim | sort >  bims
    ls *.fam | sort >  fams
    paste beds bims fams >  mergelist
    plink --merge-list mergelist  --make-bed --out raw
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
       "fill_in_bim.py ${params.output_align} $strand $manifest $inbim  raw.bim"
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


 //  Check if some of the entries are known errors to be removed at this point
 // set up delete_cmd and remove to be used next

 if (null_values.contains(mask)) {
    mask_ch = Channel.fromPath("false")
    mask_2_ch = Channel.fromPath("false")
  } else {
   mask_ch = Channel.fromPath(mask)
   mask_2_ch = Channel.fromPath(mask)
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
     file("*.{bed,bim,log,badsnps}") into aligned_ch
     file ("aligned.fam") into (fix_fam_ch, orig_fam_ch)
   script:
    base = bed.baseName
    refBase = ref.baseName
    ref_parm = ref_parm.trim()
    opt = "--a2-allele $ref $ref_parm"
    if (refBase=="empty") opt="--keep-allele-order"
    if (null_values.contains(mask)) {
	delete_cmd = ""
	remove     = " "
     } else {
	remove = " --remove mask.inds "
	delete_cmd = "list_error_inds.py ${mask_type} ${mask} '${idpat}' $fam mask.inds"
    }
    """
    ${delete_cmd}
    plink --bfile $base $opt $remove  --flip $flips --make-bed --out $output
    mv ${output}.fam aligned.fam
    grep Impossible ${output}.log | tr -d . | sed 's/.*variant//'  > ${output}.badsnps
    touch  ${output}.badsnps
    """
 }


  if (null_values.contains(params.replicates))
    replicates_ch = Channel.fromPath("False")
  else
    replicates_ch = Channel.fromPath(params.replicates)


  if (params.newpat) {
    sheet_parm = " --newpat '${params.newpat}' "
  } else {
    sheet_parm = ' '
  }



  process fixFam {
   input:
     file(samplesheet)
     file(fam) from fix_fam_ch
     file(replicates) from replicates_ch
     file(masks) from mask_2_ch
     file(sheet_columns) from sheet_cols_ch
   publishDir params.output_dir, pattern: "${output}.fam", \
             overwrite:true, mode:'copy'
   output:
     file("${output}.fam") into fixedfam_ch
   when:
     !null_values.contains(params.samplesheet)
   script:
     if (null_values.contains(params.replicates))
       replicate_parm = ""
     else 
       replicate_parm = "--replicates $replicates"
     """
     sheet2fam.py --idpat '$idpat' ${sheet_parm} --sheet-columns $sheet_columns \
                  $replicate_parm --mask $masks \
                  $samplesheet $fam $output  
     """ 
 } 


  if (null_values.contains(params.samplesheet)) {
    fixedfam_ch.combine(aligned_ch).subscribe { 
      files =  it.collect { fn -> fn.getName().replace(".*/","") }
      println "The output can be found in ${params.output_dir}: files are ${files}"
    }

 }  else {
   process publishFam {
     input:
      file (fam) from orig_fam_ch
     output:
      file("${output}.fam")
     publishDir params.output_dir, pattern: "${output}.fam", \
             overwrite:true, mode:'copy'
     script:
      " cp $fam ${output}.fam "
   }

  }
  
