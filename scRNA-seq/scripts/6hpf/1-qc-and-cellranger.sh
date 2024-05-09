working_dire="/home/fanghl/work/sdu/zhangyizhuang/rbm24/run-20240319/workflow/scRNA-seq"
trim_dire="/home/comm/sft/Trimmomatic-v0.39-d20230904"
script_dire=${working_dire}/scripts
raw_data=${working_dire}/raw-data/6hpf
res_dire=${working_dire}/6hpf/1-qc-and-cellranger
cellranger="/home/shaom/sft/cellranger-v7.2.0-d20230901"
ref_data="/home/shaom/data"


for ((i=5; i<9; i++))
do
    java -jar ${trim_dire}/trimmomatic-0.39.jar \
        SE \
        ${raw_dire}/J2309049${i}_S1_L001_R2_001.fastq.gz \
        ${res_dire}/1-filtered-data/J2309049${i}_S1_L001_R2_001.fastq.gz \
        AVGQUAL:30

    python3 ${script_dire}/utilities/align-R1-to-R2.py \
        ${raw_dire}/J2309049${i}_R1_001.fastq.gz \
        ${res_dire}/1-filtered-data/J2309049${i}_S1_L001_R2_001.fastq.gz \
        ${res_dire}/1-filtered-data/J2309049${i}_S1_L001_R1_001.fastq.gz
done


for ((i=5; i<9; i++))
do
    source ${cellranger}/sourceme.sh && cellranger count \
        --id J2309049${i} \
        --fastqs ${res_dire}/1-filtered-data \
        --sample J2309049${i} \
        --transcriptome ${cellranger_ref} \
        --output-dir ${res_dire}/2-cellranger-res \
        --localcores 30 \
        --include-introns false
done

