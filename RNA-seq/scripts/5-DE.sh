working_dire="/home/fanghl/work/sdu/zhangyizhuang/rbm24/run-20240319/workflow/RNA-seq"
script_dire=${working_dire}/script
count_dire=${working_dire}/4-count-mtx
de_dire=${working_dire}/5-DE
ref_data=${working_dire}/ref-data

time_point=("4cell" "sphere" "24HPF")

for point in "${time_point[@]}"
do
    Rscript ${script_dire}/run_DESeq2.R \
        ${count_dire}/${point}-count-mtx.tsv \
        ${count_dire}/sample.tsv \
        ${de_dire}/${point}-res.csv \
        ${de_dire}/${point}-deseq2-norm-mtx.csv

    python3 norm.py \
        --fasta_file ${ref_data}/GCF_000002035.6_GRCz11_rna_from_genomic.fna \
        --out ${de_dire}/${point} \
        ${de_dire}/${point}-count-mtx.tsv

    python3 intergrate_table.py \
        --gene_len ${de_dire}/${point}_gene_length.tsv \
        --raw_count ${de_dire}/${point}-count-mtx.tsv \
        --fpkm_count ${de_dire}/${point}_fpkm.tsv \
        --tpm_count ${de_dire}/${point}_tpm.tsv \
        --deseq_count ${de_dire}/${point}-deseq2-norm-mtx.csv \
        --deseq_res ${de_dire}/${point}-res.csv \
        --out ${point}_union.tsv
done

python3 ${script_dire}/9.2_plot.py \
    --out 4cell-sphere-24hpf \
    --norm tpm \
    4cell_union.tsv \
    sphere_union.tsv \
    24HPF_union.tsv \
    4cell sphere 24hpf
