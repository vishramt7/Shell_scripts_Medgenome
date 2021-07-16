#! /usr/bin/bash
# This script will take the report from NCBI Genome Remapping Service,
# https://www.ncbi.nlm.nih.gov/genome/tools/remap# 
# The output will be a list of markers mapped against GRCh38 (hg38)

for i in $(ls $1/*.xls)
do
    # echo $i
    # awk 'BEGIN {FS="\t"} {print $0; exit}' ${i} | sed 's/\t/\n/g' | awk '{print $1, NR}'
    # awk 'BEGIN {FS="\t"} {if ($4==$5) print $5,$13}' ${i} >> $1/report_mapped.dat # This checks the source_id and mapped_id column
    # awk 'BEGIN {FS="\t"} {print $4, $8, $5, $13}' ${i} >> $1/report_mapped.dat

done

sort -k1,1 -k2,2 -u $1/report_mapped.dat > $1/report_mapped_uniq_hg38.dat
grep -w -f $1/report_mapped_uniq_hg38.dat SARGAM_ChromNo_ChromPos.txt > report_SARGAM.intersection
