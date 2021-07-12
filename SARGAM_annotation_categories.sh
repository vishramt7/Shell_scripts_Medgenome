#! /usr/bin/bash

input_file=SARGAM_annotation_categories.tsv

#search_val="CPIC"
#echo "${search_val}"
echo "Total lines are `wc -l ${input_file}`"
awk 'BEGIN{FS="\t"} {print $0 ; exit}' ${input_file} | sed 's/\t/\n/g' | awk 'tolower($1)~/varid|description/{print $1, NR}'

awk 'BEGIN{FS="\t"} {if (NR>1) print $9,$12}' ${input_file} | awk '$2~/CPIC/ {print $1, $2}' | awk 'BEGIN{FS="_"} {if ($1=="MT") print "M",$2 ; else print $1,$2}' | sort -k1,1 -k2,2 -u > CPIC_Chrno_Chrpos.txt
# This grep helps to capture only those Chrno and Chrpos which are present in the SARGAM array.
grep -w -f CPIC_Chrno_Chrpos.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt | wc -l

awk 'BEGIN{FS="\t"} {if (NR>1) print $9,$12}' ${input_file} | awk 'tolower($2)~/blood/ {print $1, $2}' | awk 'BEGIN{FS="_"} {if ($1=="MT") print "M",$2 ; else print $1,$2}' | sort -k1,1 -k2,2 -u > Blood_Chrno_Chrpos.txt
grep -w -f Blood_Chrno_Chrpos.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt | wc -l

awk 'BEGIN{FS="\t"} {if (NR>1) print $9,$12}' ${input_file} | awk 'tolower($2)~/fingerprint/ {print $1, $2}' | awk 'BEGIN{FS="_"} {if ($1=="MT") print "M",$2 ; else print $1,$2}' | sort -k1,1 -k2,2 -u > Fingerprint_Chrno_Chrpos.txt
grep -w -f Fingerprint_Chrno_Chrpos.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt | wc -l
grep -w -f Fingerprint_Chrno_Chrpos.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt > Fingerprint_SARGAM.mapped

awk 'BEGIN{FS="\t"} {if (NR>1) print $9,$12}' ${input_file} | awk 'tolower($2)~/hla/ {print $1, $2}' | awk 'BEGIN{FS="_"} {if ($1=="MT") print "M",$2 ; else print $1,$2}' | sort -k1,1 -k2,2 -u > HLA_Chrno_Chrpos.txt
grep -w -f HLA_Chrno_Chrpos.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt | wc -l

awk 'BEGIN{FS="\t"} {if (NR>1) print $9,$12}' ${input_file} | awk 'tolower($2)~/pharm/ {print $1, $2}' | awk 'BEGIN{FS="_"} {if ($1=="MT") print "M",$2 ; else print $1,$2}' | sort -k1,1 -k2,2 -u > Pharm_Chrno_Chrpos.txt
grep -w -f Pharm_Chrno_Chrpos.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt | wc -l

awk 'BEGIN{FS="\t"} {if (NR>1) print $9,$12}' ${input_file} | awk '$2~/Mt/ {print $1, $2}' | awk 'BEGIN{FS="_"} {if ($1=="MT") print "M",$2 ; else print $1,$2}' | sort -k1,1 -k2,2 -u > Mt_Chrno_Chrpos.txt
grep -w -f Mt_Chrno_Chrpos.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt | wc -l

awk 'BEGIN{FS="\t"} {if (NR>1) print $9,$12}' ${input_file} | awk 'tolower($2)~/sars|ace/ {print $1, $2}' | awk 'BEGIN{FS="_"} {if ($1=="MT") print "M",$2 ; else print $1,$2}' | sort -k1,1 -k2,2 -u > SARS_ACE_Chrno_Chrpos.txt
grep -w -f SARS_ACE_Chrno_Chrpos.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt | wc -l

awk 'BEGIN{FS="\t"} {if (NR>1 && $6=="Blood Group") print $1}' GSA_DNAFingerprinting_BloodGroup.txt | awk 'BEGIN{FS=":"} {print $1,$2}' > GSA_blood_group.txt
awk 'BEGIN{FS="\t"} {if (NR>1 && $6=="Fingerprint SNP") print $1}' GSA_DNAFingerprinting_BloodGroup.txt | awk 'BEGIN{FS=":"} {print $1,$2}' > GSA_Fingerprint_SNP.txt

sort -k1,1 -k2,2 GSA_Fingerprint_SNP.txt | uniq > GSA_Fingerprint_SNP_uniq.txt
sort -k1,1 -k2,2 GSA_blood_group.txt | uniq > GSA_blood_group_uniq.txt

awk '{if ($1!=$2) ; else print $2,$4}' report_GSA_blood_group_uniq_mapped.txt > report_GSA_blood_group_mapped_GRCh38.txt
grep -w -f report_GSA_blood_group_mapped_GRCh38.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt | wc -l

$ awk '{if ($1!=$2) ; else print $2,$4}' report_GSA_Fingerprint_SNP_uniq_mapped.txt > report_GSA_Fingerprint_SNP_mapped_Grch38.txt
grep -w -f report_GSA_Fingerprint_SNP_mapped_Grch38.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt | wc -l
grep -w -f report_GSA_Fingerprint_SNP_mapped_Grch38.txt GWAS_catalog/SARGAM_ChromNo_ChromPos.txt > GSA_Fingerprint_SNP_SARGAM.mapped

awk 'BEGIN{FS="("}{ if ($1!="") print $1}' ACMG_genes.txt | sed 's/[[:space:]]*$//' > ACMG_gene_names.txt
grep -x -f ACMG_gene_names.txt SARGAM_gene_names.dat | sort -k1,1 | uniq | wc -l
