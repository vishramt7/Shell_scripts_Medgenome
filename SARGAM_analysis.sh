#! /usr/bin/bash

input_file=$1
no_of_columns=$(awk 'BEGIN{FS="\t"} {print NF; exit}' ${input_file})

awk '{print $0; exit}' ${input_file} | sed 's/\t/\n/g' | awk '{if ($1=="ALT") print NR}'

echo "The input file is ${input_file}, no of columns are ${no_of_columns}"

#awk '{if (NR>1) print $1,$4}' sargam.VariantClass.txt | sort -k1,1 | awk 'BEGIN{OFS="\t"} {print $1, $2}' > sargam.VariantClass_temp.dat
#awk 'BEGIN{FS="\t"} {print $3}' sargam.bestvarclass.txt | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2, $1}' > output_sargam.bestvarclass.dat
# awk '{print $0; exit}' sargam.VariMAT2.Annotation.Reduced.txt | sed 's/\t/\n/g' | awk '{if ($1=="G1000_Overall_af") print NR}'
# awk 'BEGIN{FS="\t"} {if (NR>1 && $1=="chr2") print $1,$2}' sargam.VariMAT2.Annotation.Reduced_formatted.txt | uniq > check1.dat
# awk 'BEGIN{FS="\t"} {if (NR>1 ) sum=sum+$2+$3+$4+$5} END {print sum}' SARGAM_Varclass_summary.txt
# awk 'BEGIN{FS="\t"} {print $0 ; exit}' sargam.VariMAT2.Annotation.Reduced_formatted.txt | sed 's/\t/\n/g' | awk 'tolower($1)~/^clinvar/ {print $1, NR}'
# awk 'BEGIN{FS="\t"} {if (NR>1 && $22!="NA") print $22}' SARGAM_clinical.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1) print $20}' SARGAM_clinical.txt | grep -v -i "NA" | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $21!="NA") print $21}' SARGAM_clinical.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $23!="NA") print $23}' SARGAM_clinical.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $24!="NA") print $24}' SARGAM_clinical.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $25!="NA") print $25}' SARGAM_clinical.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $26!="NA" && $26!="NA:NA") print $26}' SARGAM_clinical.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $27!="NA") print $27}' SARGAM_clinical.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $28!="NA") print $28}' SARGAM_clinical.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $29!="NA") print $29}' SARGAM_clinical.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $3!="NA" && $3!=".") print $3}' SARGAM_ancestry.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $2!="NA") print $2}' SARGAM_COSMIC_REFSEQ_GWAS.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $2!="NA") print $2}' SARGAM_Allele_frequency.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $3!="NA") print $3}' SARGAM_COSMIC_REFSEQ_GWAS.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $4!="NA") print $4}' SARGAM_COSMIC_REFSEQ_GWAS.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $2!="NA") print $2}' SARGAM_COSMIC_REFSEQ_GWAS.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $2!="NA") print $2}' SARGAM_SNP.txt | wc -l
# awk 'BEGIN{FS="\t"} {if (NR>1 && $8!="NA") print $8 }' SARGAM_Allele_frequency.txt | awk '{if ($1<0.01) sum++} END{print sum}'
# awk 'BEGIN{FS="\t"} {if (NR>1 && $8!="NA") print $8 }' SARGAM_Allele_frequency.txt | awk '{if ($1>0.05) sum++} END{print sum}'
# awk 'BEGIN{FS="\t"} {if (NR>1 && $8!="NA") print $8 }' SARGAM_Allele_frequency.txt | awk '{if ($1>=0.01 && $1<=0.05) sum++} END{print sum}'
# awk 'BEGIN{FS="\t"} {if (NR>1 && $10!="NA") print $10 }' SARGAM_Allele_frequency.txt | awk '{if ($1<0.01) sum++} END{print sum}'
# awk 'BEGIN{FS="\t"} {if (NR>1 && $10!="NA") print $10 }' SARGAM_Allele_frequency.txt | awk '{if ($1>0.05) sum++} END{print sum}'
# awk 'BEGIN{FS="\t"} {if (NR>1 && $10!="NA") print $10 }' SARGAM_Allele_frequency.txt | awk '{if ($1>=0.01 && $1<=0.05) sum++} END{print sum}'


