#!/usr/bin/bash 


for i in $(cut -f 1 ../positions.txt | sort -nk1,1 | uniq)
do
	echo ${i}
	/BXRX_DS_ANALYSIS5/plink2 --bgen /BXRX_DS_RESOURCES/DATA_FROM_UKBIOBANK/ARRAY_DATASET/IMPUTED/ukb_imp_chr${i}_v3.bgen ref-first --sample /BXRX_DS_RESOURCES/DATA_FROM_UKBIOBANK/ARRAY_DATASET/IMPUTED/ukb42406_imp_chr${i}_v3_s487320.sample --extract range ../positions.txt --recode vcf --out UKBB_imput_chr${i}
	#ls /BXRX_DS_RESOURCES/DATA_FROM_UKBIOBANK/ARRAY_DATASET/IMPUTED/ukb_imp_chr${i}_v3.bgen /BXRX_DS_RESOURCES/DATA_FROM_UKBIOBANK/ARRAY_DATASET/IMPUTED/ukb42406_imp_chr${i}_v3_s487320.sample positions_chr${i}.txt UKBB_imput_chr${i}

done
