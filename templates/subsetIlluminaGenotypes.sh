#
#   Generate a random subset of an Illumina GWAS data set (gsgt format)
#   ===================================================================
#
#   This script creates a random subset of a cohort genotype dataset in
#   the illumina genotype genome studio (gtgs) format. Alongside the
#   standard illumina genotype data files:
#	+ sample report (csv)
#	+ locus report (csv)
#	+ genotype reports (csv.gz)
#   it takes as input the desired number of snvs and desired number of
#   samples and outputs a fully functional subset that can be converted
#   into plink format (for example) and carried through a typical GWAS
#   pipeine.
#
#   The script creates a subset of the samples by randomly sampling
#   records from the illumina sample report. It creates a subset of
#   the snvs by randomly sampling from the illumina locus report.
#   Finally, it then uses these to subset the genotype reports so that
#   the output data set is fully consistent and operational with plink.
#
#
#######################################################################

#
#   Input arguments
#

number_of_genotype_reports=10
number_of_samples_in_subset=100
number_of_snvs_in_subset=1000
input_data_prefix='WTS_H3Africa_Wonkam_2018.04'
output_data_prefix=''

#
#   Workflow definition
#

function workflow {

    custom_call prepare_workspace || { exit 1; }

    custom_call get_random_subset_of_sample_report || { exit 1; }

    custom_call get_random_subset_of_locus_report || { exit 1; }

    for (( i=1; i<=$number_of_genotype_reports; i++ ))
    do
        custom_call get_matching_subset_of_genotype_report "$i" || { exit 1; }
    done
}


#
#   Process definitions
#

function custom_call {
    local method=$1
    shift

    printf "running *$method* with arguments ("$@")..."
    if $method "$@" ; then
        echo "...done"
    else
        echo "...failed!"
        return 1
    fi
}

function prepare_workspace {
    mkdir -p data-subset
    mkdir -p temp
}


function get_random_subset_of_sample_report {

    local sample_report=${input_data_prefix}_DNA_Summary_Report.csv
    local sample_report_of_subset=data-subset/${output_data_prefix}sample-report.csv

    create_file $sample_report_of_subset || { return 1; }

    copy_file_header $sample_report $sample_report_of_subset || { return 1; }

    copy_random_subset_of_records \
        $sample_report \
        $sample_report_of_subset \
        $number_of_samples_in_subset \
        || { return 1; }

    get_list_of_ids $sample_report_of_subset data-subset/sample-ids.txt || { return 1; }
}

function get_random_subset_of_locus_report {

    local locus_report=${input_data_prefix}_Locus_Summary_Report.csv
    local locus_report_of_subset=data-subset/${output_data_prefix}locus-report.csv

    create_file $locus_report_of_subset || { return 1; }

    copy_file_header $locus_report $locus_report_of_subset || { return 1; }

    copy_random_subset_of_records \
        $locus_report \
        $locus_report_of_subset \
        $number_of_snvs_in_subset \
        || { return 1; }

   get_list_of_ids $locus_report_of_subset data-subset/locus-ids.txt || { return 1; } 
}

function get_matching_subset_of_genotype_report {
    local genotype_report_id=$1

    local genotype_report=${input_data_prefix}_gtReport_File-${genotype_report_id}.csv.gz
    local genotype_report_of_subset=data-subset/${output_data_prefix}genotype-report-${genotype_report_id}.csv
    local genotype_report_header

    create_file $genotype_report_of_subset || { return 1; }

    get_file_header $genotype_report genotype_report_header || { return 1; }

    modify_file_header \
        "$genotype_report_header" \
        genotype_report_header_modified \
        "s/Content/Content (random subset: ${number_of_samples_in_subset} samples x ${number_of_snvs_in_subset} snvs)/g" \
        || { return 1; }

    write_file_header \
        "$genotype_report_header_modified" \
        $genotype_report_of_subset \
        || { return 1; }

    filter_by_sample \
        $genotype_report \
        data-subset/sample-ids.txt \
        temp/${genotype_report_id} \
        || { return 1; }

    filter_by_snv \
        temp/${genotype_report_id} \
        data-subset/locus-ids.txt \
        $genotype_report_of_subset \
        || { return 1; }

    gzip -f $genotype_report_of_subset

    clean_temporary_files  || { return 1; }
}

function create_file {
    local new_file=$1

    touch $new_file
}

function copy_file_header {
    local source_file=$1
    local target_file=$2

    head -n 1 $source_file >> $target_file
}

function copy_random_subset_of_records {
    local source_file=$1
    local target_file=$2
    local number_of_records_to_copy=$3

    awk '(NR>1) {print}' $source_file \
        | sort -R \
        | head -n $number_of_records_to_copy \
        >> $target_file
}

function get_list_of_ids {
    local report_file=$1
    local output_list=$2

    awk -F ',' '{print $2}' $report_file > $output_list
}

function get_file_header {
    local source_file=$1
    local target_buffer=$2

    local number_of_lines_in_genotype_report_header=10

    eval $target_buffer="'$(zcat $source_file | head -n $number_of_lines_in_genotype_report_header)'"
}

function modify_file_header {
    local source_buffer=$1
    local target_buffer=$2
    local replacement_pattern=$3

    eval $target_buffer="'$(echo "$source_buffer" | sed "${replacement_pattern}")'"
}

function write_file_header {
    local source_buffer=$1
    local target_file=$2

    echo "$source_buffer" >> $target_file
}

function filter_by_sample {
    local source_report=$1
    local list_of_samples_to_keep=$2
    local target_report=$3

    zcat $source_report \
        | awk -F ',' 'FNR==NR {a[$1];next} $2 in a' \
        $list_of_samples_to_keep \
        - \
        > $target_report
}

function filter_by_snv {
    local source_report=$1
    local list_of_snvs_to_keep=$2
    local target_report=$3

    awk -F ',' 'FNR==NR {a[$1];next} $1 in a' \
        $list_of_snvs_to_keep \
        $source_report \
        >> $target_report
}

function clean_temporary_files {

    rm temp/${genotype_report_id}
}


#
#   Run the workflow
#

workflow
