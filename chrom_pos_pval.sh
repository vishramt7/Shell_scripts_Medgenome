#!/usr/bin/bash
# This script will extract the chrom pos and pval for an input study 

#infile=$1

#outfile=$2

#pvalue_column=$(awk 'BEGIN{FS="\t"} {for (i=1; i<=NF; i++) { if (tolower($i) ~ /all_inv_var_meta_p|p_value|pval/) print i}; exit}' ${infile} | head -n1)
#awk -v col=$pvalue_column 'BEGIN {OFS=FS="\t"}NR>1{ print $1,$2,$col}' ${infile} > ${outfile}.txt

for i in $(ls *.txt)
do
        file=$(basename ${i} .txt)
        echo ${file}
	#wc -l ${i}
        #pvalue_column=$(awk '{for (i=1; i<=NF; i++) { if (tolower($i) ~ /all_inv_var_het_p/) print i}; exit}' ${i} | head -n1)
        #awk '{if ($3<0.00000005) print $1,$2}' ${i} > ${file}_gws.chrmo_pos
        #awk -v col=$pvalue_column 'NR>1 {if ($col<0.00001) print $1,$2}' ${i} > ${file}_sugests.chrmo_pos
        #sort -nk1,1 -nk2,2 ${i}        | uniq >> non_UKBB_sugests.chrmo_pos
        #./rsids_extract_vcf_gz.pl ${i}
	#sort -nk1,1 -nk2,2 ${file}_gws.chrmo_pos | uniq >> concat_gws.chrmo_pos
done

#sort -nk1,1 -nk2,2 concat_gws.chrmo_pos | uniq -c | sort -nrk1,1 | awk '{if ($1==6) print $2,$3}' > common_6.chr_pos

