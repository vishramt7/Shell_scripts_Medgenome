#!/usr/bin/bash 

for i in $(ls *b37.txt)
do 
	file=$(basename ${i} .txt)
	echo ${file}
	Rscript manhattan_modified.R ${file}
done
