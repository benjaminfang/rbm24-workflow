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

countmtx \
    --count_level gene_name\
    --out ${count_dire}/4cell-count-mtx \
    ${read_anno}/M-4cell-*-assigned-* \
    ${read_anno}/Sibling-4cell-*-assigned-*

countmtx \
    --count_level gene_name\
    --out ${count_dire}/sphere-count-mtx \
    ${read_anno}/M-sphere-*-assigned-* \
    ${read_anno}/Sibling-sphere-*-assigned-*

countmtx \
    --count_level gene_name\
    --out ${count_dire}/24HPF-count-mtx \
    ${read_anno}/M-24HPF-*-assigned-* \
    ${read_anno}/Sibling-24HPF-*-assigned-*

