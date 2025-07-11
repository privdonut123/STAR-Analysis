#!/bin/bash

input_dir=$(ls -d input/track_qa/*/)


for dir in $input_dir
do
    rapidity_substr=$(echo $dir | cut -d'/' -f3 | cut -c5- | tr '_' '.')
    
    # echo $rapidity_substr
    root4star -l -q fitProjections.C\(\"$dir\"\)
    root4star -l -q pt_qa.C\(\"$rapidity_substr\"\)
    mv output/pt_qa.root output/res_qa_eta_$rapidity_substr.root
done