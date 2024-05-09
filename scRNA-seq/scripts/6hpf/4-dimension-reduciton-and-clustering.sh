working_dire="/home/fanghl/work/sdu/zhangyizhuang/rbm24/run-20240319/workflow/scRNA-seq"
script_dire=${working_dire}/scripts
data_dire=${working_dire}/6hpf
res_dire=${working_dire}/6hpf/4-dimension-reduction-and-barcodes-clustering

python3 ${script_dire}/utilities/dimension-reduciton.py \
    --germ_marker ${script_dire}/germcell-markers.txt \
    --mt_genes ${script_dire}/mt-genes.txt \
    --max_scale_value 10 \
    --out ${res_dire}/s95_s96_s97_s98 \
    ${data_dire}/2-filter-of-barcodes/2-integrated-data/s95_s96_s97_s98 \
    ${data_dire}/3-selection-of-features/s95_s96_s97_s98_hvg.txt

python3 ${script_dire}/utilities/clustering.py \
    --out ${res_dire}/s95_s96_s97_s98 \
    --revolution .6 \
    ${res_dire}/s95_s96_s97_s98.h5ad
