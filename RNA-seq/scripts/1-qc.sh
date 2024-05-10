working_dire="/home/fanghl/work/git-repo/rbm24-workflow/RNA-seq"
raw_data=${working_dire}/raw-data
clean_data=${working_dire}/1-clean-data

for item in "$raw_data"/*
do
    if [[ $item == *"L1_1"* ]]
    then
        echo "-------------------"
        item_clean=$(basename $item)
        item_clean=${clean_data}/"${item_clean/L1_1/L1_1-clean}"
        item2="${item/L1_1/L1_2}"
        item2_clean=$(basename $item2)
        item2_clean=${clean_data}/"${item2_clean/L1_2/L1_2-clean}"
        echo $item
        echo $item_clean
        echo $item2
        echo $item2_clean

        cutadapt \
            -a AGATCGGAAGAGCACACG \
            -A AGATCGGAAGAGCGTCGT \
            --minimum-length 50 \
            -q 15,10 \
            --trim-n \
            -o ${item_clean} \
            -p ${item2_clean} \
            $item \
            $item2
    fi
done 

