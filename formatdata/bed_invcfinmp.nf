
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

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths;
import sun.nio.fs.UnixPath;
import java.security.MessageDigest;


// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}



def helps = [ 'help' : 'help' ]
allowed_params = ['file_bed', 'file_vcf', 'file_ref']

//params.file_vcf=""
params.input_dir=""
params.input_pat=""

inpat = "${params.input_dir}/${params.input_pat}"

bedinitial=Channel.create()
biminitial=Channel.create()
Channel
   .fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName }  \
      .ifEmpty { error "No matching plink files" }        \
      .map { a -> [checker(a[1]), checker(a[2]), checker(a[3])] }\
      .separate(bedinitial, biminiral) { a -> [a,a[2]] }

#plink --bfile /spaces/wenlongchen/version3/qc_dnachip_round02/output/gwas_input --geno 0.01 --hwe 0.0001 --snps-only --maf 0.01 --mind 0.01 --make-bed --out temp1 --threads 12

process extractrsname{
  input :
    file(bim) from 
  

}
process convertInVcf {
   memory plink_mem_req
   cpus max_plink_cores
   input :
    file(bed), file(bim), file(fam) from bedbeforefilt
   output :
    base=bed.baseName
    file("${base}.vcf")  into vcfi
   script:
     base= plink[0].baseName
     """
     """
 }

//plink --bfile ${base} --threads ${max_plink_cores} --recode vcf --out $base --keep-allele-order --snps-only --threads $max_plink_cores
//plink --bfile temp2 --exclude multiallelic_snps_remove.txt --make-bed --out temp3 --threads 12

