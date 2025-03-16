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

process extract_site_from_vcf {
    tag "extract_site_${target_name}_${site_name}"
    label "bigmem"

    input:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), val(site_name), file(site_file)
    output:
        tuple val(target_name), val(chrm), val(start), val(end), file(tag_target_vcf), val(site_name), file(site_file)
    script:
        base = "${target_name}_${site_name}"
        tag_target_vcf = "${base}.bcf"
        """
        tabix -f ${target_vcf}
        bcftools view --regions-file ${site_file} ${target_vcf}  -Ob -o ${tag_target_vcf}
        """
}

process sort_vcf {
    tag "sort_${target_name}"
    label "bigmem"

    input:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), val(site_name), file(site_file)
    output:
        tuple val(target_name), val(chrm), val(start), val(end), file(sorted_vcf), val(site_name), file(site_file)
    script:
        base = file(target_vcf.baseName).baseName
        sorted_vcf = "${base}_sorted.bcf"
        """
        bcftools sort ${target_vcf} -T . -Ob -o ${sorted_vcf}
        """
}

process fill_tags_VCF {
    tag "fill_tags_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf), val(chrm)
    output:
        tuple val(dataset), file(out_vcf), val(chrm)
    script:
        base = file(vcf.baseName).baseName
        out_vcf = "${base}_AF.vcf.gz"
        """
        tabix -f ${vcf}
        bcftools +fill-tags ${vcf}  -Oz -o ${out_vcf}
        """
}

process vcf_freq {
    tag "freq_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf), val(chrm)
    output:
        tuple val(dataset), file(frq), val(chrm)
    script:
        base = file(vcf.baseName).baseName
        frq = "${base}.frq"
        """
        tabix -f ${vcf}
        echo -e 'CHROM\tPOS\tSNP\tREF\tALT\tAF' > ${frq}
        bcftools annotate  --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' ${vcf} | \
        bcftools query  -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/AF\\n' >> ${frq}
        """
}

process vcf_map {
    tag "map_${target_name}"
    label "bigmem"

    input:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), val(site_name), file(site_file)
    output:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), file("${target_vcf}.csi"), val(site_name), file(map)
    script:
        base = file(target_vcf.baseName).baseName
        map = "${base}.map"
        """
        tabix -f ${target_vcf}
        bcftools annotate  --set-id +'%CHROM\\_%POS' ${target_vcf} | \
        bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' > ${map}
        """
}

process vcf_map_simple {
    tag "map_${target_name}"
    label "bigmem"

    input:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), val(site_name), file(site_file)
    output:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), file("${target_vcf}.csi"), val(site_name), file(map)
    script:
        base = file(target_vcf.baseName).baseName
        map = "${base}.map"
        """
        tabix ${target_vcf}
        bcftools query -f '%CHROM\\t%POS\\n' ${target_vcf} > ${map}
        """
}

process info_freq {
    tag "info_freq_${dataset}_${ref_panel}"
    label "bigmem"

    input:
        tuple val(dataset), file(info), val(ref_panel)
    output:
        tuple val(dataset), val(ref_panel), file(frq), file(info)
    script:
        base = info.baseName
        frq = "${base}.frq"
        """
        awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6}' ${info} > ${frq}
        """
}

process get_vcf_sites {
    tag "get_sites_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${tagName}"
    label "bigmem"

    input:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(target_vcf), val(tagName)
    output:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(sites), val(tagName)
    script:
        base = file(target_vcf.baseName).baseName
        sites = "${base}.sites"
        """
        tabix -f ${target_vcf}
        echo -e 'SNP' > ${sites}
        bcftools annotate  --set-id +'%CHROM\\_%POS' ${target_vcf} | \
        bcftools query  -f '%ID\\n' >> ${sites}
        """
}

process generate_chunks {
    tag "generate_chunks_${target_name}"
    label "small"

    input:
        tuple val(target_name), file(mapFile), val(chunk_size)
    output:
        tuple val(target_name), file(chunkFile)
    script:
        if(params.imp_chunk){ chunk = params.imp_chunk } else{ chunk='' }
        chromosomes = ''
        chunkFile = "chunks.txt"
        template "generate_chunks.py"
}

process generate_chunks_vcf {
    tag "generate_chunks_${target_name}"
    label 'imputation'                                                          
    memory params.high_memory
    input:
        tuple val(target_name), file(vcf), file(mapFile), val(chrms), val(chunk_size)
    output:
        tuple val(target_name), file(vcf), file(chunkFile)
    script:
        if(params.imp_chunk){ chunk = params.imp_chunk } else{ chunk='' }
        chromosomes = chrms
        chunkFile = "chunks.txt"
        template "generate_chunks.py"
}


process split_target_to_chunk {
    tag "split_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${tagName}"
    label 'imputation'                                                          
    memory params.high_memory   
    //label "bigmem"
    maxForks 30

    input:
        tuple val(target_name), val(chrm), val(chunk_start), val(chunk_end), val(tagName), file(tag_target_vcf)
    output:
        tuple val(target_name), val(chrm), val(chunk_start), val(chunk_end), val(tagName), file(vcf_chunk_out)
    script:
        base = file(tag_target_vcf.baseName).baseName
        vcf_chunk_out = "${base}.${chrm}_${chunk_start}-${chunk_end}_${tagName}.bcf"
        start = chunk_start - params.imp_buffer_size
        if(chunk_start.toInteger() - params.imp_buffer_size.toInteger() <= 0){ end = 1 }
        end = chunk_end.toInteger() + params.imp_buffer_size.toInteger()
        """
        tabix ${tag_target_vcf}
        bcftools view --regions ${chrm}:${start}-${end} ${tag_target_vcf} -Ob -o ${vcf_chunk_out}
        """
}


process extract_pop {
    tag "extract_pop_vcf_${target_name}_${chrm}_${ref_name}"
    // publishDir "${params.outDir}/imputed/vcfs/${ref_name}/${prefix}/${target_name}", overwrite: true, mode:'link', pattern: '*vcf.gz*'
    label "bigmem"

    input:
        tuple val(chrm), val(target_name), val(ref_name), val(tagName), file(vcf), file(vcf_tbi), file(sample_to_extract), val(prefix)
    output:
        tuple val(chrm), val(target_name), val(ref_name), val(tagName), file(vcf_out), file("${vcf_out}.tbi")
    script:
        base = file(vcf.baseName).baseName
        vcf_out = "${base}_${prefix}.vcf.gz"
        """
        bcftools \
            view \
             \
            --force-samples \
            --samples-file ${sample_to_extract} \
            ${vcf} \
            -Oz -o ${vcf_out} 
        tabix ${vcf_out} 
        """
}
