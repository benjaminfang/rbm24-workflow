working_dire="/home/fanghl/work/git-repo/rbm24-workflow/RNA-seq"
mapped_data=${working_dire}/2-mapped-data
ref_data=${working_dire}/ref-data/
read_anno=${working_dire}/3-reads-anno

for item in "$mapped_data"/*
do
    if [[ $item == *"fixmate-sort-redup-sort.sam" ]]
    then
        echo "-------------------"
        basename=$(basename $item)
        rdpos=${read_anno}/${basename/-fixmate-sort-redup-sort.sam/}
        echo $item
        echo $basename
        echo $rdpos

        annoread.py \
            --fasta_file ${ref_data}/GCF_000002035.6_GRCz11_rna_from_genomic.fna \
            --gtf_file ${ref_data}/GCF_000002035.6_GRCz11_genomic.gtf \
            --read_type pe \
            --out ${rdpos} \
            ${item}
    fi
done 

