working_dire="/home/fanghl/work/git-repo/rbm24-workflow/RNA-seq"
clean_data=${working_dire}/1-clean-data
mapped_data=${working_dire}/2-mapped-data
ref_data=${working_dire}/ref-data


for item in "$clean_data"/*
do
    if [[ $item == *"L1_1"* ]]
    then
        echo "-------------------"
        item2="${item//L1_1/L1_2}"
        echo $item
        echo $item2
        item_basename=${mapped_data}/$(basename $item)
        map_res="${item_basename/_L1_1-clean.fq.gz/.sam}"
        map_res_fixmate="${item_basename/_L1_1-clean.fq.gz/-fixmate.bam}"
        map_res_fixmate_sort="${item_basename/_L1_1-clean.fq.gz/-fixmate-sort.bam}"
        map_res_fixmate_sort_redup="${item_basename/_L1_1-clean.fq.gz/-fixmate-sort-redup.bam}"
        map_res_fixmate_sort_redup_sort="${item_basename/_L1_1-clean.fq.gz/-fixmate-sort-redup-sort.sam}"
        echo $map_res
        echo $map_res_fixmate
        echo $map_res_fixmate_sort
        echo $map_res_fixmate_sort_redup
        echo $map_res_fixmate_sort_redup_sort

        bowtie2 \
            -x ${ref_data}/rna-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic.fna \
            -1 ${item} \
            -2 ${item2} \
            -S ${map_res} \
            -k 3 \
            --threads 8

        samtools fixmate \
            --threads 8 \
            -m \
            -O bam \
            ${map_res} \
            ${map_res_fixmate}

        samtools sort \
            --threads 8 \
            -O bam \
            -o ${map_res_fixmate_sort} \
            ${map_res_fixmate}

        samtools markdup \
            -O bam \
            --threads 8 \
            -r ${map_res_fixmate_sort} \
            ${map_res_fixmate_sort_redup}

        samtools sort \
            --threads 8 \
            -n \
            -O sam \
            -o ${map_res_fixmate_sort_redup_sort} \
            ${map_res_fixmate_sort_redup}
    fi
done 

