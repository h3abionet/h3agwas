#!/usr/bin/env nextflow

/*
========================================================================================
=                                 h3achipimputation                                    =
========================================================================================
 h3achipimputation imputation functions and processes.
----------------------------------------------------------------------------------------
 @Authors

----------------------------------------------------------------------------------------
 @Homepage / @Documentation
  https://github.com/h3abionet/chipimputation
----------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------
*/


// check if files exist [name, file1, file2, ...]
def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
}

def check_chromosome_vcf(dataset, dataset_vcf, map_file, chrms){
    """
    Check if specified chromosomes exist in VCF file
    """
    chromosomes_ = [:]
    chromosomes_['ALL'] = []
    valid_chrms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
    not_chrs = []
    in_chrs = []
    notValid_chrs = []
    chromosomes_[dataset] = file(map_file).readLines().collect{ it.split('\t')[0] }//.unique().sort()
    chromosomes_[dataset].each { chrm ->
        chrm = chrm.trim()
        // if(!(chrm in chromosomes_['ALL'])) { //todo not working
            if (chrm in valid_chrms){
                chromosomes_['ALL'] << chrm
            }
            else{
                notValid_chrs << chrm
            }
        // }
    }
    if (chrms == '' || chrms == 'ALL'){
        chromosomes = chromosomes_['ALL']
    }
    else{
        chrms.toString().split(',').each { chrm ->
            if (!(chrm in chromosomes_['ALL'])){
                not_chrs << chrm
            }
            else{
                in_chrs << chrm
            }
        }
        chromosomes = in_chrs
        if (in_chrs.isEmpty()){
            exit 1,  "|-- ERROR- No Chromosome(s) found not in target(s) dataset(s)! The pipeline will exit."
        }
        if (!(not_chrs.isEmpty())){
            log.info "|-- WARN- Chromosome(s) ${not_chrs.join(', ')} not in target datasets and will be ignored."
        }
    }
    return [ dataset, file(dataset_vcf), map_file, in_chrs ]
}


def get_chromosome_vcf(vcf){
    """
    Check if specified chromosomes exist in VCF file
    """

    chrom = file(vcf).readLines().unique().collect{ it as int }.sort()
    
    return in_chrs
}


/*
 * Check user's provided chromosomes vs those in map file
 */
process check_chromosome {
    tag "check_chromosome_${target}"
    input:
        tuple val(target), file(target_vcf)
    output:
        tuple val(target), file(chrom_file)
    script:
        base = file(target_vcf.baseName).baseName
        chrom_file = "${base}_chromosomes.txt"
        """
        zcat ${target_vcf} | grep -v "^#" | awk -F' ' '\$1 ~ /^[0-9]+\$/ {print \$1}' | sort -n | uniq >  ${chrom_file}
        """
}

process get_chromosome {
    tag "get_chromosome_${dataset}"
    label 'imputation'                                                          
    memory params.high_memory 

    input:
        tuple val(dataset), file(dataset_vcf)
    output:
        tuple val(dataset), file(dataset_vcf), file(chrom_file)
    script:
        base = file(dataset_vcf.baseName).baseName
        chrom_file = "${base}_chromosomes.txt"
        """
        bcftools query -f '%CHROM\\t%POS\\n' ${dataset_vcf} > ${chrom_file}
        """
}

process split_vcf_chromosome {
    tag "split_vcf_chrm_${dataset}"
    label 'imputation'                                                          
    memory params.high_memory 
    input:
        tuple val(dataset), file(dataset_vcf)
    output:
        tuple val(dataset), file("${base}*.vcf.gz")
    script:
        base = file(dataset_vcf.baseName).baseName
        """
        bcftools view ${dataset_vcf} -Ou -o ${dataset_vcf.baseName}
        SnpSift split ${dataset_vcf.baseName}
        N=1
        for FILE in ${base}.*.vcf
        do
            bcftools view \${FILE} -Oz -o ${base}.\$N.vcf.gz
            rm \${FILE}
            N=\$((N+1))
        done
        """
}

process split_vcf_chunk {
    tag "split_vcf_${dataset}_${chunk_size}"
    label 'imputation'                                                          
    memory params.high_memory   
    input:
        tuple val(dataset), file(dataset_vcf), val(chunk_size)
    output:
        tuple val(dataset), file("${base}*.vcf.gz")
    script:
        base = file(dataset_vcf.baseName).baseName
        """
        SnpSift split -l ${chunk_size} ${dataset_vcf}
        N=1
        for FILE in ${base}.*.vcf
        do
            bcftools view \${FILE} -Oz -o ${base}.\$N.vcf.gz
            rm \${FILE}
            N=\$((N+1))
        done
        """
}

process check_mismatch {
    tag "check_mismatch_${target}_${chrm}_${start}_${end}"
    //label "medium"
    label 'imputation'
    
    input:
        tuple val(target), val(chrm), val(start), val(end), file(target_vcf), file(reference_genome)
    output:
        tuple val(target), val(chrm), val(start), val(end), file(target_vcf), file("${base}_checkRef_warn.log"), file("${base}_checkRef_summary.log")
    script:
        base = file(target_vcf.baseName).baseName
        """
        bcftools norm  --check-ref w -f ${reference_genome} ${target_vcf} -Oz -o /dev/null
        cp .command.err ${base}_checkRef_warn.log
        bcftools +fixref ${target_vcf} -- -f ${reference_genome} 2>&1 | tee "${base}_checkRef_summary.log"
        rm -f ${base}_clean_mind.*
        """
}

def no_mismatch(target_name, warn, summary){
    mismatch = 0
    file(warn).readLines().each{ it ->
        if(it.contains("REF_MISMATCH")){
            mismatch += 1
        }
    }
    if ( mismatch != 0 ) {
        System.err.println "|-- ERROR: ${mismatch} ref mismatch sites found in '${target_name}' dataset! The pipeline will exit. Check this file ${summary}."
        exit 1
    }
    return mismatch
}

process target_qc {
    tag "target_qc_${target_name}"
    label "bigmem"

    input:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf)
    output:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_qc_vcf)
    script:
        base = file(target_vcf.baseName).baseName
        target_qc_vcf = "${base}_qc.bcf"
        """
        bcftools norm --rm-dup both ${base}_noALT.bcf -Ob -o ${target_qc_vcf}
        rm -f ${base}_noALT.snp ${base}_noALT.bcf
        """
}

process qc_dupl {
    tag "dupl_qc_${dataset}_${chrm}_${start}_${end}"
    label 'imputation'
    memory params.high_memory

    input:
        tuple val(dataset), val(chrm), val(start), val(end), file(dataset_vcf)
    output:
        tuple val(dataset),val(chrm), val(start), val(end), file(dataset_qc_vcf)
    script:
        base = file(dataset_vcf.baseName).baseName
        dataset_qc_vcf = "${base}_${chrm}_${start}_${end}_qc.bcf"
        """
        tabix -f ${dataset_vcf}
        bcftools norm  --rm-dup both ${dataset_vcf} -Ob -o ${dataset_qc_vcf}
        tabix -f ${dataset_qc_vcf}
        """
}

process split_multi_allelic {
    tag "split_multi_${dataset}_${chrm}_${start}_${end}"
    label 'imputation'                                                          
    memory params.high_memory                                                   
                                                  

    input:
        tuple val(dataset), val(chrm), val(start), val(end), file(vcf)
    output:
        tuple val(dataset), val(chrm), val(start), val(end), file(vcf_out)
    script:
        base = file(vcf.baseName).baseName
        vcf_out = "${base}.norm.bcf"
        """
        bcftools norm -m- ${vcf} -Ob -o ${vcf_out}
        tabix -f ${vcf_out}
        """
}

process fill_tags_vcf {
    tag "fill_tags_${dataset}_${chrm}"
    label 'imputation'                                                          
    memory params.high_memory                                                   
    input:
        tuple val(dataset), val(chrm), val(start), val(end), file(vcf)
    output:
        tuple val(dataset), val(chrm), val(start), val(end), file(vcf_out)
    script:
        vcf_out = "${vcf.baseName}_af.bcf"
        """
        bcftools +fill-tags ${vcf}  -Ob -o ${vcf_out}
        tabix -f ${vcf_out}
        """
}

process filter_min_ac {
    tag "min_ac_${dataset}_${chrm}_${start}_${end}"
    label 'imputation'                                                          
    memory params.high_memory  
    input:
        tuple val(dataset), val(chrm), val(start), val(end), file(vcf), val(params)

    output:
        tuple val(dataset), val(chrm), val(start), val(end), file(vcf_out), file(map_file)

    script:
        base = file(vcf.baseName).baseName
        vcf_out = "${base}_ac.bcf"
        map_file = "${base}_ac.chrom.txt"
        """
        bcftools view ${params} ${vcf} -Ob -o ${vcf_out}
        tabix -f ${vcf_out}
        bcftools query -f '%CHROM\\t%POS\\n' ${vcf_out} > ${map_file}
        """
 }


process qc_site_missingness {
    tag "site_missingness_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}_${tagName}"
    label "bigmem"

    input:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(impute_vcf), file(imputed_info), val(tagName), val(qc_params)
    output:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(impute_vcf_qc), file(imputed_info), val(tagName)
    script:
        base = file(impute_vcf.baseName).baseName
        impute_vcf_qc = "${base}_qc.vcf.gz"
        """
        vcftools --gzvcf ${impute_vcf} --keep-INFO-all ${qc_params} --recode --stdout | bgzip > ${impute_vcf_qc}
        tabix ${impute_vcf_qc}
        """
}

process sites_only {
    tag "sites_only_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}_${tagName}"
    label "bigmem"

    input:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(impute_vcf_qc), file(imputed_info), val(tagName)

    output:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(sites_vcf), file(imputed_info), val(tagName)

    script:
        base = base = file(impute_vcf_qc.baseName).baseName
        sites_vcf = "${base}_sites.vcf.gz"
        """
        tabix ${impute_vcf_qc}
        bcftools view ${impute_vcf_qc} --drop-genotypes  -Oz -o ${sites_vcf}
        tabix ${sites_vcf}
        """
}

process combine_vcfs_chrm {
   tag "combine_${chrm}_${target_name}_${ref_name}_${tagName}"
   publishDir "${params.output_dir}/imputed/vcfs/${ref_name}/all/${target_name}/${tagName}", overwrite: true, mode:'copy', pattern: '*vcf.gz*'
   label "bigmem"
   
   input:
       tuple val(chrm), val(target_name), val(ref_name), val(tagName), val(vcfs)
   output:
       tuple val(chrm), val(target_name), val(ref_name), val(tagName), file(vcf_out), file("${vcf_out}.tbi")
   script:
       vcf_out = "${target_name}_${ref_name}_${tagName}_chr${chrm}.vcf.gz"
       if(vcfs.size() > 1){
            """
            bcftools concat ${vcfs.join(' ')} --allow-overlaps  | \
            bcftools sort -T .  -Oz -o ${vcf_out}
            tabix ${vcf_out}
            """
       }
       else if(vcfs.size() == 1){
            """
            cp ${vcfs.join(' ')} ${vcf_out}
            tabix ${vcf_out}
            """
       }
}


process combine_vcfs {
   tag "combine_${dataset}_${ref_name}"
   label "bigmem"
   
   input:
       tuple val(dataset), val(ref_name), val(vcfs)
   output:
       tuple val(dataset), val(ref_name), file(vcf_out)
   script:
       vcf_out = "${dataset}_${ref_name}_phased.vcf.gz"
       if(vcfs.size() > 1){
            """
            bcftools concat ${vcfs.join(' ')} --allow-overlaps  | \
            bcftools sort -T .  -Oz -o ${vcf_out}
            tabix ${vcf_out}
            """
       }
       else if(vcfs.size() == 1){
            """
            cp ${vcfs.join(' ')} ${vcf_out}
            tabix ${vcf_out}
            """
       }
}

// """
// Index VCF
// """
// process index_vcf {
//     tag "index_${dataset}_${chrm}"
//     publishDir "${params.output_dir}/imputed/vcfs/${ref_name}/${prefix}/${target_name}", overwrite: true, mode:'copy', pattern: '*vcf.gz.tbi'
//     label "medium"
   
//     input:
//         tuple dataset, file(vcf), chrm

//     output:
//         tuple dataset, file(vcf), file(index)

//     script:
//         index = "${vcf}.tbi"
//         """
//         tabix ${vcf}
//         """
// }

"""
Combine impute info chunks to chromosomes
"""
process combine_infos {
    tag "combine_infos_${target_name}_${ref_name}_${tagName}"
    publishDir "${params.output_dir}/imputed/infos/${ref_name}", overwrite: true, mode:'copy', pattern: '*imputed_info'
    label "medium"

    input:
        tuple val(target_name), val(ref_name), val(infos), val(tagName)
    output:
        tuple val(target_name), val(ref_name), file(comb_info), val(tagName)
    script:
        comb_info = "${target_name}_${ref_name}.info"
        """
        head -n1 ${infos[0]} > ${comb_info}
        tail -q -n +2 ${infos.join(' ')} >> ${comb_info}
        """
}


"""
Combine csvs
"""
process combine_csvs {
    tag "combine_csvs_${target_name}_${ref_name}"
    publishDir "${params.output_dir}/imputed/${ref_name}", overwrite: true, mode:'copy', pattern: "*${out_ext}"
    label "medium"

    input:
        tuple val(target_name), val(ref_name), val(csvs), val(out_ext)
    output:
        tuple val(target_name), val(ref_name), file(comb)
    script:
        comb = "${target_name}_${ref_name}.${out_ext}"
        """
        head -n1 ${csvs[0]} > ${comb}
        tail -q -n +2 ${csvs.join(' ')} >> ${comb}
        """
}

// def helpMessage() {
//     log.info"""
//     =========================================
//     h3achipimputation v${params.version}
//     =========================================
//     Usage:

//     The typical command for running the pipeline is as follows:

//     nextflow run h3abionet/chipimputation --reads '*_R{1,2}.fastq.gz' -profile standard,docker

//     Mandatory arguments (Must be specified in the configuration file, and must be surrounded with quotes):
//       --target_datasets             Path to input study data (Can be one ou multiple for multiple runs)
//       --genome                      Human reference genome for checking REF mismatch
//       --ref_panels                  Reference panels to impute to (Can be one ou multiple for multiple runs)
//       -profile                      Configuration profile to use. Can use multiple (comma separated)
//                                     Available: standard, conda, docker, singularity, test

//     Other options:
//       --outDir                      The output directory where the results will be saved
//       --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
//       --name                        Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
//       --project_name                Project name. If not specified, target file name will be used as project name
//     """.stripIndent()
// }
