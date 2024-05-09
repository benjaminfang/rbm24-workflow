working_dire="/home/fanghl/work/sdu/zhangyizhuang/scRNA-seq/run-20240319/workflow/scRNA-seq"
cellranger_dire="/home/comm/sft/cellranger-7.2.0"
ref_data=${working_dire}/ref-data

source ${cellreanger_dire}/sourceme.bash && cellranger mkref \
    --genome=${ref_data}/cellranger-zebrafish-ref \
    --fasta=${ref_data}/GCF_000002035.6_GRCz11_genomic.fna \
    --genes=${ref_data}/GCF_000002035.6_GRCz11_genomic.gtf

