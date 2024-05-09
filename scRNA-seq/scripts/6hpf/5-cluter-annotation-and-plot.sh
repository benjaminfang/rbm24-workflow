working_dire="/home/fanghl/work/sdu/zhangyizhuang/rbm24/run-20240319/workflow/scRNA-seq"
script_dire=${working_dire}/scripts
data_dire=${working_dire}/6hpf/4-dimension-reduction-and-barcodes-clustering
res_dire=${working_dire}/6hpf/5-annotation-of-clusters


python3 ${script_dire}/utilities/plot_res.py \
    --ctl s97,s98 \
    --exp s95,s96 \
    --marker_genes ${res_dire}/cluster-marker-genes-6hpf \
    --germ_marker ${script_dire}/germcell-markers.txt \
    --count_dire ${working_dire}/6hpf/2-filter-of-barcodes/2-integrated-data/s95_s96_s97_s98_nm \
    --heatmap_genes ${res_dire}/heatmap_genes \
    --scatter_genes ${res_dire}/scatter_genes \
    --violin_genes ${res_dire}/violin_genes \
    --out ${res_dire}/6hpf_figures \
    ${data_dire}/s95_s96_s97_s98_rv_0.2_cell_cluster.tsv \
    ${res_dire}/cluster-info-6hpf 

