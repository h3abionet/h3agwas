#
#   Generate a random subset of an
#   Illumina GWAS data set in gsgt
#
######################################################################

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

    head -n 1 $sample_report > $sample_report_of_subset
    awk '(NR>1) {print}' $sample_report \
        | sort -R \
        | head -n $number_of_samples_in_subset \
        >> $sample_report_of_subset
    awk -F ',' '{print $2}' $sample_report_of_subset \
        >> data-subset/sample-ids.txt
}

function get_random_subset_of_locus_report {
    local locus_report=${input_data_prefix}_Locus_Summary_Report.csv
    local locus_report_of_subset=data-subset/${output_data_prefix}locus-report.csv

    head -n 1 $locus_report > $locus_report_of_subset
    awk '(NR>1) {print}' $locus_report \
       | sort -R \
        | head -n $number_of_snvs_in_subset \
        >> $locus_report_of_subset
    awk -F ',' '{print $2}' $locus_report_of_subset \
        >> data-subset/locus-ids.txt
}

function get_matching_subset_of_genotype_report {
    local genotype_report_id=$1
    local genotype_report=${input_data_prefix}_gtReport_File-${genotype_report_id}.csv.gz
    local genotype_report_of_subset=data-subset/${output_data_prefix}genotype-report-${genotype_report_id}.csv
    local genotype_report_header

    genotype_report_header=$(zcat $genotype_report | head -n 10)
    echo "$genotype_report_header" \
        | sed "s/Content/Content (random subset: ${number_of_samples_in_subset} samples x ${number_of_snvs_in_subset} snvs)/g" \
        > $genotype_report_of_subset
    zcat $genotype_report \
        | awk -F ',' 'FNR==NR {a[$1];next} $2 in a' data-subset/sample-ids.txt - \
        > temp/${genotype_report_id}
    awk -F ',' 'FNR==NR {a[$1];next} $1 in a' data-subset/locus-ids.txt temp/${genotype_report_id} \
        >> $genotype_report_of_subset
    gzip -f $genotype_report_of_subset
    rm temp/${genotype_report_id}
}

#
#   Run the workflow
#

workflow
