working_dire="/home/fanghl/work/sdu/zhangyizhuang/rbm24/run-20240319/workflow/scRNA-seq/"
script_dire=${working_dire}/scripts
data_dire=${working_dire}/6hpf/1-qc-and-cellranger/2-cellranger-res
res_dire=${working_dire}/6hpf/2-filter-of-barcodes

for ((i=5; i<9; i++))
do
    python3 ${script_dire}/utilities/filter-barcodes.py \
        --mt_file ${script_dire}/mt-genes.txt \
        --decontx_script ${script_dire}/utilities/decontX.R \
        --out ${res_dire}/1-barcodes-filtered/s9${i} \
        ${data_dire}/J2309049${i}
done

python3 ${script_dire}/utilities/intergrate-datasets.py \
    --mt_genes ${script_dire}/mt-genes.txt \
    --germ_marker ${script_dire}/germcell-markers.txt \
    --out ${res_dire}/2-integrated-data \
    ${res_dire}/1-barcodes-filtered/s95 \
    ${res_dire}/1-barcodes-filtered/s96 \
    ${res_dire}/1-barcodes-filtered/s97 \
    ${res_dire}/1-barcodes-filtered/s98 \
    ${data_dire}/J23090495 \
    ${data_dire}/J23090496 \
    ${data_dire}/J23090497 \
    ${data_dire}/J23090498

python3 ${script_dire}/utilities/normalize_data_by_1e4.py \
    --out_dire ${res_dire}/2-integrated-data/s95_s96_s97_s98_nm \
    ${res_dire}/2-integrated-data/s95_s96_s97_s98

python3 ${script_dire}/utilities/stats-and-plot-barcode-info.py \
    --stats_file ${res_dire}/2-integrated-data/s95_s96_s97_s98_barcodes_stat.txt \
    --out ${res_dire}/density-cumulation-6hpf \
    s95 s96 s97 s98
