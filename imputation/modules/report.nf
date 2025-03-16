#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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


// process quality_control
//     tag "perform QC on the imputed files"

process filter_info {
    tag "filter_${dataset_name}_${tagName}_${ref_panels.join('-')}"
    label "bigmem"
    label 'imputation'                                                          

    input:
        tuple val(dataset_name), val(ref_panels), val(ref_infos), val(tagName), val(info_cutoff)
    output:
        tuple val(dataset_name), val(ref_panels), file("${well_out}_${info_cutoff}.tsv"), file("${acc_out}_${info_cutoff}.tsv"), val(tagName)
    script:
        comb_info = "${dataset_name}_${tagName}_${ref_panels.join('-')}.imputed_info"
        well_out = "${comb_info}_well_imputed"
        acc_out = "${comb_info}_accuracy"
        infos = ref_infos.join(',')
        impute_info_cutoff = params.impute_info_cutoff
        template "filter_info_minimac.py"
}


process report_site_by_maf {
    tag "site_by_maf_${dataset_name}"
    label "bigmem"
    label 'imputation'                                                          

    input:
        tuple val(dataset_name), file(sites)
    output:
        tuple val(dataset_name), file("${site_by_maf}_number_of_imputed.tsv"), file("${site_by_maf}_number_of_imputed_summary.tsv"), file("${site_by_maf}_avg_info.tsv")
    script:
        site_by_maf = "${sites.baseName}_report_by_maf"
        group = 'REF_PANEL'
        template "report_well_imputed.py"
}

"""
Report 1: Generate a file of well imputed snps by maf for a dataset (key) for all reference panels (datasets/refpanels)
"""
//TODO do this by chromosomes for each dataset
process report_well_imputed_by_target {
    tag "report_wellImputed_${target_name}_${ref_panels.split(',').join('-')}"
    publishDir "${params.output_dir}/imputed/reports/${ref_panels}", overwrite: true, mode:'copy'
    label 'imputation'                                                          
    label "medium"
    
    input:
        tuple val(target_name), val(ref_panels), file(inWell_imputed)
    
    output:
        tuple val(target_name), val(ref_panels), file("${out_prefix}.tsv"), file("${out_prefix}_summary.tsv")
    
    script:
        out_prefix = "${inWell_imputed.baseName}.imputed_info_performance_by_maf_report"
        // outWell_imputed = "${target_name}_${ref_panels}_${chrms}.imputed_info_performance_by_maf_report"
        group = "REF_PANEL"
        template "report_well_imputed.py"
}

"""
Plot performance all reference panels by maf for a dataset
"""
process plot_performance_target{
    tag "plot_performance_dataset_${target_name}_${ref_panels}_${chrms}"
    publishDir "${params.output_dir}/imputed/plots/${ref_panels}", overwrite: true, mode:'copy'
    label 'imputation'                                                          
    
    input:
        tuple val(target_name), val(ref_panels), file(well_imputed_report), file(well_imputed_report_summary), val(group)
    output:
        tuple val(target_name), val(ref_panels), file(plot_by_maf)
    script:
        plot_by_maf = "${well_imputed_report.baseName}.pdf"
        report = well_imputed_report
        xlab = "MAF bins"
        ylab = "Number of well imputed SNPs"
        template "plot_results_by_maf.R"
}

"""
Repor 2: Accuracy all reference panels by maf for a dataset
"""
process report_accuracy_target {
    tag "report_acc_${target_name}_${ref_panels}"
    publishDir "${params.output_dir}/imputed/reports/${ref_panels}/", overwrite: true, mode:'copy'
    label 'imputation'                                                          
    label "medium"
    input:
        tuple val(target_name), val(ref_panels), file(inSNP_acc), val(group)
    output:
        tuple val(target_name), val(ref_panels), file(outSNP_acc), val(group)
    script:
        outSNP_acc = "${inSNP_acc.baseName}.imputed_info_report_accuracy.tsv"
        template "report_accuracy_by_maf.py"
}

"""
Plot accuracy all reference panels by maf for a dataset
"""
process plot_accuracy_target{
    tag "plot_accuracy_dataset_${target_name}_${ref_panels}"
    publishDir "${params.output_dir}/imputed/plots/${ref_panels}", overwrite: true, mode:'copy'
    label 'imputation'                                                          
    input:
        tuple val(target_name), val(ref_panels), file(accuracy_report), val(group)
    output:
        tuple val(target_name), val(ref_panels), file(plot_by_maf)
    script:
        plot_by_maf = "${accuracy_report.baseName}_accuracy_by_maf.pdf"
        report = accuracy_report
        xlab = "MAF bins"
        ylab = "Concordance rate"
        template "plot_results_by_maf.R"
}

"""
Step: generate allele frequency
"""
process generate_frequency {
    tag "frq_${target_name}_${ref_name}"
    publishDir "${params.output_dir}/imputed/frqs/${ref_name}", overwrite: true, mode:'copy', pattern: '*frq'
    label "medium"
    label 'imputation'                                                          
    input:
        tuple val(target_name), val(ref_name), file(impute_vcf),  file(ref_vcf)
    output:
        tuple val(target_name), val(ref_name), file(dataset_frq), file(ref_frq)
    script:
        dataset_frq = "${file(impute_vcf.baseName).baseName}.frq"
        ref_frq = "${file(ref_vcf.baseName).baseName}.frq"
        
        """
        # For datastet
        echo -e 'CHR\tPOS\tSNP\tREF\tALT\tAF' > ${dataset_frq}
        bcftools view -m2 -M2 -v snps ${impute_vcf} | bcftools query -f '%CHROM\t%POS\t%CHROM\\_%POS\\_%REF\\_%ALT\t%REF\t%ALT\t%INFO/AF\\n' >> ${dataset_frq}
        
        # For the reference panel
        echo -e 'CHR\tPOS\tSNP\tREF\tALT\tAF' > ${ref_frq}
        bcftools view -m2 -M2 -v snps ${ref_vcf} | bcftools +fill-tags -Oz -o ${ref_name}_AF.vcf.gz -- -t AF
        bcftools query -f '%CHROM\t%POS\t%CHROM\\_%POS\\_%REF\\_%ALT\t%REF\t%ALT\t%INFO/AF\\n' ${ref_name}_AF.vcf.gz >> ${ref_frq}
        
        """
}

"""
Plot number of imputed SNPs over the mean r2 for all reference panels
"""
process plot_r2_SNPpos {
    tag "plot_r2_SNPpos_${target_name}_${ref_name}_${chrm}"
    publishDir "${params.output_dir}/imputed/plots/${ref_name}", overwrite: true, mode:'copy'
    label 'imputation'                                                          
    label "medium"
    input:
    tuple val(target_name), val(ref_name), file(target_info), val(maf_thresh), file(target_frq)
    output:
    tuple val(target_name), val(ref_name), file(output)
    script:
    info = target_info
    target = target_frq
    output = "${target_name}_${ref_name}_r2_SNPpos.pdf"
    template "r2_pos_plot.R"
}


"""
Plot frequency of imputed SNPs against SNP frequencies in reference panels
"""
process plot_freq_comparison {
    tag "plot_freq_comparison_${target_name}_${ref_name}"
    publishDir "${params.output_dir}/imputed/plots/${ref_name}/freq_comparison", overwrite: true, mode:'copy'
    label 'imputation'                                                          
    label "medium"
    input:
        tuple val(target_name), val(ref_name), file(target_info), file(target_frq), file(ref_frq) 
    output:
        tuple val(target_name), val(ref_name), file(outputcolor)
    script:
        info = target_info
        target = target_frq
        frq = ref_frq
        //output = "${target_name}_${ref_name}_${chrm}_freq_comparison.pdf"
        outputcolor = "${target_name}_${ref_name}_freq_comparison_color.pdf"
        template "AF_comparison.R"
}
       

"""
Plot number of imputed SNPs over the mean r2 for all reference panels
"""
process plot_r2_SNPcount {
    tag "plot_r2_SNPcount_${target_name}_${ref_panels}"
    publishDir "${params.output_dir}/imputed/plots/${ref_panels}", overwrite: true, mode:'copy'
    label 'imputation'                                                          
    label "medium"
    input:
        tuple val(target_name), val(ref_panels), val(infos)
    output:
        tuple val(target_name), val(ref_panels), file(plot_out) 
    script:
    plot_out = "${target_name}_${ref_panels}_r2_SNPcount.pdf"
    impute_info_cutoff = params.impute_info_cutoff
    template "r2_Frequency_plot.R"
}

"""
Plot histograms of number of imputed SNPs over the mean r2 for all reference panels
"""
process plot_hist_r2_SNPcount {
    tag "plot_hist_r2_SNPcount_${target_name}_${ref_panels}"
    publishDir "${params.output_dir}/imputed/plots/${ref_panels}/", overwrite: true, mode:'copy'
    label 'imputation'                                                          
    input:
        tuple val(target_name), val(ref_panels), val(infos)
    output:
        tuple val(target_name), val(ref_panels), file(plot_out)
    script:
    plot_out = "${target_name}_${ref_panels}_r2_SNPcount_hist.pdf"
    impute_info_cutoff = params.impute_info_cutoff
    template "r2_Frequency_plot_histogram.R"
}

"""
Plot MAF of imputed SNPs over r2 for all references
"""
process plot_MAF_r2 {
    tag "plot_MAF_r2_${target_name}_${ref_panels}"
    publishDir "${params.output_dir}/imputed/plots/${ref_panels}", overwrite: true, mode:'copy'
    label "medium"
    label 'imputation'                                                          
    input:
        tuple val(target_name), val(ref_panels), val(infos) 
    output:
        tuple val(target_name), val(ref_panels), file(plot_out) 
    script:
    plot_out = "${target_name}_${ref_panels}_MAF_r2.pdf"
    impute_info_cutoff = params.impute_info_cutoff
    template "Frequency_r2_MAF_plot.R"
}

process average_r2 {
    tag "average_r2_${target_name}_${ref_panels}"
    publishDir "${params.output_dir}/imputed/rsquared/${ref_panels}", overwrite: true, mode:'copy'
    label 'imputation'                                                          
    input:
        tuple val(target_name), val(ref_panels), val(file_infos) 
    output:
        tuple val(target_name), val(ref_panels), file(meanr2_out) 
    script:
    meanr2_out = "${target_name}_${ref_panels}"
    template "average_r2.R"
}
