#!/usr/bin/bash 

for i in $(ls *UKBB*)
do 
	file=$(basename ${i} .txt)
	echo ${file}
	Rscript manhattan_modified.R ${file}
	cp ${file}_manhattan.png ../Manhattan_plots/
done
