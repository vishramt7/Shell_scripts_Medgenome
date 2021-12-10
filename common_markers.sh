#!/usr/bin/bash 


# Path of genotyped ukbb data
ukbb_geno_path=/BXRX_DS_RESOURCES/DATA_FROM_UKBIOBANK/ARRAY_DATASET/GENOTYPES
# Path of MAC filtered 1000G data
G1000_path=/BXRX_ANALYSIS_DATA/ANALYSIS4/GWAS_POPLTN_STUDY/GSA/GA100k_1000G_plan1/beagle/data
for i in {2..22}
do 
	echo ${i}
        plink --bfile ${ukbb_geno_path}/ukb_cal_chr${i}_v2 --recode vcf --out UKBB_geno_${i}
	
	vcftools --vcf UKBB_geno_${i}.vcf --gzdiff ${G1000_path}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_mac_hwe.vcf.recode.vcf.gz --diff-site --out ukbb_v_1000g_chr${i}

	awk 'BEGIN{OFS="\t"}{if ($4=="B") print $1,$2,$2,"id"}' ukbb_v_1000g_chr${i}.diff.sites_in_files > positions_chr${i}.txt

	plink --vcf UKBB_geno_${i}.vcf --extract range positions_chr${i}.txt --recode vcf --out ukbb_common_chr${i}
	plink --vcf ${G1000_path}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_mac_hwe.vcf.recode.vcf.gz --extract range positions_chr${i}.txt --recode vcf --out 1000g_common_chr${i}
	rm UKBB_geno_${i}.vcf
	gzip ukbb_common_chr${i}.vcf 1000g_common_chr${i}.vcf
done
