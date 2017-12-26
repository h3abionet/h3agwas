/* (C) H3ABionet
 GPL

 Scott Hazelhrust
*/

params.vcf_merge_table="/dataB/aux/RsMergeArch.bcp.gz"
params.reference      ="/dataB/aux/Homo_sapiens.GRCh37.75.dna.toplevel.fa"
params.dbsnp_all      ="/dataB/aux/37/All_20170403.vcf.gz"
params.manifest       ="/dataB/aux/HumanOmni5-4-v1-0-D.csv"
params.chipdescription = "/external/diskA/novartis/metadata/HumanOmni5-4v1_B_Physical-and-Genetic-Coordinates.txt"
inpat = "${params.input_dir}/${params.input_pat}"
params.output         = "chip"
params.strandreport    = false

params.vcf = false

plink_src = Channel.create()

if (params.strandreport) {
  strand_ch = Channel.fromPath(params.strandreport);
} else {
  File empty = new File("empty.txt")
  empty.createNewFile()
  strand_ch = Channel.fromPath("empty.txt")
}
	 



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
     for (c in channels[1..-1]){
         start = start.phase(c) .map(combineElts)
     }
     return start;
}




  output      = params.output


  rsmerge_ch  = Channel.create()

  samplesheet = Channel.fromPath(params.samplesheet)
  ref_ch      = Channel.fromPath(params.reference)
  remap       = Channel.fromPath("scripts/remap_map_and_ped.pl")
  dbsnp_all   = Channel.fromPath(params.dbsnp_all)
  array        = Channel.fromPath(params.chipdescription)
  array2        = Channel.fromPath(params.chipdescription)
  manifest_csv   = Channel.fromPath(params.manifest)
  report        = Channel.fromPath(inpat)

  chroms = 1..26



  

  process illumina2lgen {
    maxForks 32
    input:
       set file(report), file(array) from report.combine(array)
    output:
       set file("${output}.ped"), file("${output}.map") into ped_ch
    script:
	abbrev = 1
        output = report.baseName
        template "topbot2plink.py"
  }




  process bedfiy {
    
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



 process combine_bed {
   input:
     file(bed) from bed_ch.toList()
     file(bim) from bim_ch.toList()
     file(fam) from fam_ch.toList()
   output:
     set file("raw.{bed,bim,fam,log}") into plink_src
     file ('orig.fam') into fix_fam_ch
   script:
    """
    ls *.bed | sort > beds
    ls *.bim | sort >  bims
    ls *.fam | sort >  fams
    paste beds bims fams >  mergelist
    plink --merge-list mergelist --make-bed --out $output
    mv ${output}.fam orig.fam
    """
 } 


 process getFlips {
   input:
     file(strandreport) from strand_ch
   output:
     file(flips) into flip_ch
   script:
    flips = "flips.txt"
    template "strandmismatch.py"
 }


 process alignStrand {
   input:
    set file(bed), file(bim), file(fam), file(log) from plink_src
    file(flips) from flip_ch
    publishDir params.output_dir, pattern: "*.{bed,bim,log}", \
        overwrite:true, mode:'copy'
   output:
    set file("$output.{bed,bim,log}") into aligned_ch
   script:
    base = bed.baseName
    """
    plink --bfile $base --flip $flips --make-bed --out $output
    cp $log combine.log
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
    abbreviate = "1"
    template "sheet2fam.py"
 } 




