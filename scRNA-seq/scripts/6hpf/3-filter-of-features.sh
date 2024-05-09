working_dire="/home/fanghl/work/sdu/zhangyizhuang/rbm24/run-20240319/workflow/scRNA-seq"
script_dire=${working_dire}/scripts
data_dire=${working_dire}/6hpf/2-filter-of-barcodes/2-integrated-data
res_dire=${working_dire}/6hpf/3-selection-of-features

python3 ${script_dire}/utilities/select-features.py \
    --germ_marker ${script_dire}/germcell-markers.txt \
    --mt_genes ${script_dire}/mt-genes.txt \
    --cell_cycle ${script_dire}/cell-cycle-house-keeping-genes.txt \
    --out ${res_dire}/s95_s96_s97_s98 \
    ${data_dire}/s95_s96_s97_s98

