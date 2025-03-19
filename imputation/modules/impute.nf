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

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

process impute_minimac4 {
    tag "imp_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}_${tagName}"
    label 'imputation'                                                          
    memory params.high_memory   
    publishDir "${params.output_dir}/imputed/${ref_name}/${target_name}", overwrite: true, mode:'copy'

    input:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), file(target_phased_vcf), file(target_phased_vcf_tbi), val(ref_name), file(ref_vcf), file(ref_m3vcf), val(tagName)
    output:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file("${base}_imputed.dose.vcf.gz"), file("${base}_imputed.info"), val(tagName), emit : out
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file("${base}_imputed.dose.vcf.gz"), file("${base}_imputed.info"), val(tagName), file(ref_vcf), emit : outwithref
    shell:
        base = "${file(target_phased_vcf.baseName).baseName}_${tagName}_${chrm}_${chunk_start}-${chunk_end}"
        """
        minimac4 \
            --refHaps ${ref_m3vcf} \
            --haps ${target_phased_vcf} \
            --format GT,DS \
            --allTypedSites \
            --minRatio ${params.imp_minRatio} \
            --chr ${chrm} --start ${chunk_start} --end ${chunk_end} --window ${params.imp_buffer_size} \
            --prefix ${base}_imputed \
            --cpus ${task.cpus}
        """
}

process impute5 {
    label 'imputation'                                                          

    tag "imp_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}_${tagName}"
    publishDir "${params.output_dir}/imputed/impute5/${ref_name}", overwrite: true, mode:'copy', pattern: '*vcf.gz'
    label "bigmem_impute5"
    input:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), file(target_phased_vcf), file(target_phased_vcf_tbi), val(ref_name), file(ref_vcf), file(ref_imp5), file(ref_imp5_idx), val(tagName), file(impute5_genetic_map)
    output:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(impute5_out), val(tagName)
    shell:
        impute5_out = "${file(target_phased_vcf.baseName).baseName}_${tagName}_${chrm}_${chunk_start}-${chunk_end}.vcf.gz"
        """
        impute5_1.1.5_static \
            --h ${ref_imp5} \
            --m ${impute5_genetic_map} \
            --g ${target_phased_vcf} \
            --r ${chrm}:${chunk_start}-${chunk_end} \
            --o ${impute5_out}
        """
}

process generate_impute5_info {
    label 'imputation'                                                          
    tag "imp_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}"
    publishDir "${params.output_dir}/imputed/impute5/${ref_name}", overwrite: true, mode:'copy'
    input:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(impute5_out)
    output:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file("${generate_info}.txt")
    shell:
        generate_info = "${target_name}_${chrm}_${chunk_start}-${chunk_end}"
        """
        bcftools annotate --set-id '%CHROM\\:%POS\\:%REF\\:%FIRST_ALT' ${impute5_out} > ${generate_info}.vcf.gz
        bcftools query -f '%ID\\t%REF\\t%FIRST_ALT\\t%INFO/AF\\t%INFO/AF\\t%INFO/INFO\\t%INFO/IMP\\n' ${generate_info}.vcf.gz > noheader_${generate_info}_draft.txt
        (echo -e "SNP\\tREF(0)\\tALT(1)\\tALT_Frq\\tMAF\\tRsq\\tGenotyped"; cat noheader_${generate_info}_draft.txt ) > ${generate_info}_test.txt
        awk 'BEGIN{ OFS="\\t" } { if (\$7=="."){ \$7="Genotyped" } print }' ${generate_info}_test.txt > ${generate_info}_replaced.txt
        awk 'BEGIN{ OFS="\\t" } { if (\$7==1){ \$7="Imputed" } print }' ${generate_info}_replaced.txt > ${generate_info}.txt
        """
}


process impute_minimac4_1 {
    tag "imp_${target_name1}_${chrm1}:${chunk_start1}-${chunk_end1}_${ref_name1}_${tagName1}"
    label 'imputation'                                                          
    memory params.high_memory   
    input:
        tuple val(chrm1), val(chunk_start1), val(chunk_end1), val(target_name1), file(target_phased_vcf1), val(ref_name1), file(ref_vcf1), file(ref_m3vcf1), val(tagName1)
    output:
        tuple val(chrm1), val(chunk_start1), val(chunk_end1), val(target_name1), val(ref_name1), file("${base}_imputed.dose.vcf.gz"), file("${base}_imputed.info"), val(tagName1)
    shell:
        base = "${file(target_phased_vcf1.baseName).baseName}"
        """
        nblines=\$(zcat ${target_phased_vcf1} | grep -v "^#" | wc -l)
        if (( \$nblines > 1 ))
        then
            bcftools annotate  --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' ${target_phased_vcf1} | 
            vcftools --vcf - --keep-INFO-all --max-missing 0.05 --hwe 0.00001 --mac 1 --recode --stdout | bgzip > ${base}.id.vcf.gz
            minimac4 \
                --cpus ${task.cpus} \
                --refHaps ${ref_m3vcf1} \
                --haps ${base}.id.vcf.gz \
                --format GT,DS \
                --allTypedSites \
                --minRatio ${params.imp_minRatio} \
                --chr ${chrm1} --start ${chunk_start1} --end ${chunk_end1} --window ${params.imp_buffer_size} \
                --prefix ${base}_imputed
        else
             touch ${base}_imputed.dose.vcf && bgzip ${base}_imputed.dose.vcf
             touch ${base}_imputed.info
        fi
        """
}

"""
Combine impute chunks to chromosomes
"""
process combineImpute {
    tag "impComb_${target_name}_${ref_name}_${chrm}"
    publishDir "${params.output_dir}/imputed/${ref_name}/${target_name}", overwrite: true, mode:'copy', pattern: '*imputed.vcf.gz'
    label "bigmem"
    label 'imputation'                                                          
    
    input:
        tuple val(target_name), val(ref_name), val(chrm), val(imputed_files)
    output:
        tuple val(target_name), val(ref_name), val(chrm), file(comb_impute)// into combineImpute
    script:
        comb_impute = "${target_name}_${ref_name}_chr${chrm}.imputed.vcf.gz"
        """
        bcftools concat ${imputed_files.join(' ')} | \
        bcftools +fill-tags -- -t AC,AN,AF,MAF | \
        bcftools sort -T . -Oz -o ${comb_impute}
        """
}

"""
Combine impute info chunks to chromosomes
"""
process combineInfo {
    tag "infoComb_${target_name}_${ref_name}_${chrm}"
    publishDir "${params.output_dir}/imputed/${ref_name}/${target_name}", overwrite: true, mode:'copy', pattern: '*imputed_info'
    label "medium"
    label 'imputation'                                                          

    input:
        tuple val(target_name), val(ref_name), val(chrm), file(info_files) //from infoCombine.values()
    output:
        tuple val(target_name), val(ref_name), val(chrm), file(comb_info) //into combineInfo_frq
    script:
        comb_info = "${target_name}_${ref_name}_chr${chrm}.imputed_info"
        """
        head -n1 ${info_files[0]} > ${comb_info}
        tail -q -n +2 ${info_files.join(' ')} >> ${comb_info}
        """
}

"""
Filtering all reference panels by maf for a dataset
"""
//TODO generate filtered info by reference panels.
process filter_info_by_target {
    tag "filter_${target_name}_${ref_panels}_${chrms}"
    publishDir "${params.output_dir}/imputed/reports/${ref_panels}", overwrite: true, mode:'copy', pattern: "${comb_info}*"
    label "medium"
    label 'imputation'
    
    input:
        tuple val(target_name), val(ref_panels), val(chrms), val(ref_infos) //from target_infos.values()
    
    output:
        tuple val(target_name), val(ref_panels), file("${out_prefix}_well_imputed.tsv"), file("${out_prefix}_accuracy.tsv")  //into target_info_Well
        // set target_name, ref_panels, file("${acc_out}.tsv") //into target_info_Acc,target_info_Acc_1
    
    script:
        chrms = chrms.split(',')[0]+"-"+chrms.split(',')[-1]
        out_prefix = "${target_name}_${ref_panels.split(',').join('-')}_${chrms}.imputed_info"
        ref_infos = ref_infos
        datasets = ref_panels
        impute_info_cutoff = params.impute_info_cutoff
        template "filter_info_minimac.py"
}

