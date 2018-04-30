#!/bin/bash

declare -a dataSets=("TT" "QCD")

for d in "${dataSets[@]}"
do
    echo "----------------------------------"
    echo "        Combining $d Trees        "
    echo "----------------------------------"
    rm output-files/${d}/make_training_trees_${d}.root
    find output-files/${d} -type f -name "*" > ${d}_temp.txt
    root -q -b -l 'CombineTree.C("'${d}'_temp.txt",  "output-files/'${d}'/make_training_trees_'${d}'.root")'
    rm ${d}_temp.txt
done
