working_dire="/home/fanghl/work/git-repo/rbm24-workflow/RNA-seq"
read_anno=${working_dire}/3-reads-anno
count_dire=${working_dire}/4-count-mtx


for item in "$read_anno"/*
do
    if [[ $item == *".rdpos" ]]
    then
        echo "-------------------"
        echo $item
        assigned=${item/.rdpos/-assigned.rdpos}
        echo $assigned

        assignread \
            --gene_name drop \
            --out $assigned \
            $item
    fi
done 

time_point=("4cell" "sphere" "24HPF")

for point in "${time_point[@]}"
do
    countmtx \
        --count_level gene_name\
        --out ${count_dire}/${point}-count-mtx \
        ${read_anno}/${point}-*-assigned-* \
        ${read_anno}/Sibling-${point}-*-assigned-*
done

