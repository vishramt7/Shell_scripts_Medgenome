#!/bin/bash

# 1. Parse arguments
# a. Parse flag, if it exists
savemethod="none" # start by assuming user doesn't want to save intermediate files
verbose=0 # default
skipibd=0 # default

while [[ $# -gt 0 ]]
do
key="$1"
case "$1" in
	--help | -h)
		echo "OPTIONS:"
		echo "-h, --help                show flag descriptions"
		echo "-l, --save-logs           saves logs from QC procedure to ./[FILENAME]_QC_logs"
		echo "-a, --save-all            saves all intermediate files from QC procedure, including logs, to ./[FILENAME]_QC_intermediate_files"
		echo "-v, --verbose             unsupresses XWAS output"
		echo "-g, --debug               saves logs, intermediate files, and keep XWAS output"
		echo "-s, --skip-ibd            skips IBD analysis step in QC procedure, should only be used if you are confident you have no relatedness in your sample or have used another method for IBD analysis"
		exit 0
		;;
	--save-logs | -l)
		savemethod="logs"
		;;
	--save-all | -a)
		savemethod="all"
		;;
	--verbose | -v)
		verbose=1
		;;
	--debug | -g)
		savemethod="all"
		verbose=1
		;;
	--skip-ibd | -s)
		skipibd=1
		;;
	*)
		;;
esac
shift
done
paramfile="$key"

# b. Check that correct arguments are given
if [ ! -f "$paramfile" ]; then
	echo "Incorrect parameter usage."
	echo "File not found:" "$paramfile"
	echo "Correct useage: ./run_QC.sh [ options ] /path/to/param_file.txt"
	echo "Type ./run_QC.sh --help for futher details."
	exit 1
fi

# c. Parse param file info
fname=$( awk '$1=="filename" {print $2}' $paramfile ) # file name without extension
xwasloc=$( awk '$1=="xwasloc" {print $2}' $paramfile ) # location of bin
binaryplink=$( awk '$1=="plinkformat" {print $2}' $paramfile ) # ped = ped/map format, bed = binary plink format

# d. If user provides ped/map format, convert to binary plink format
if [ $binaryplink = "ped" ]; then
	${xwasloc}/xwas --file ${fname} --make-bed --out ${fname} --silent
fi

# Removing indels, by extracting only snps
plink --noweb --bfile ${fname} --make-bed --snps-only --out ${fname}_snp_only		# This is the only step which uses plink executable !!
# Removing rs enteries separated by ;. Taking only one of the rsids
awk 'BEGIN{OFS="\t"; FS="[\t;]"} {print $1,$(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF}' ${fname}_snp_only.bim > temp.bim
rm ${fname}_snp_only.bim
mv temp.bim ${fname}_snp_only.bim

# 2. Perform pre-QC
perl ${xwasloc}/get_independent_individials.pl ${fname}_snp_only

if [ "$verbose" -eq 1 ]; then
	${xwasloc}/xwas --noweb --bfile ${fname}_snp_only --remove ${fname}_snp_only.PC_Sib_clean.ind_to_remove --make-bed --out ${fname}.PC_sib_clean
else
	${xwasloc}/xwas --noweb --bfile ${fname}_snp_only --remove ${fname}_snp_only.PC_Sib_clean.ind_to_remove --make-bed --out ${fname}.PC_sib_clean --silent
fi

perl ${xwasloc}/set_parental_information_to_be_missing.pl ${fname}.PC_sib_clean.fam ${fname}.preprocessed.fam
mv ${fname}.PC_sib_clean.bed ${fname}.preprocessed.bed
mv ${fname}.PC_sib_clean.bim ${fname}.preprocessed.bim

# 3. Create temporary parameter file
touch temp_params.txt
echo "filename ${fname}.preprocessed" > temp_params.txt
sed '1d' ${paramfile} >> temp_params.txt
old_fname="$fname" # save name to delete files later

# 4. Run QC
paramfile="temp_params.txt"
if [ "$verbose" -eq 1 ]; then
	xwasQ="xwas"
else
	xwasQ="xwas --silent"
fi

echo -e "3: Parsing parameter file"
fname=$( awk '$1=="filename" {print $2}' $paramfile )		# name of dataset file (without the extension)
xwasloc=$( awk '$1=="xwasloc" {print $2}' $paramfile )		# location of extended plink executable
eigstratloc=$( awk '$1=="eigstratloc" {print $2}' $paramfile ) # smartpca and convertf file locations.
exclind=$( awk '$1=="exclind" {print $2}' $paramfile )		# 0 if no file exists, 1 if file exists
xchrpca=$( awk '$1=="excludexchrPCA" {print $2}' $paramfile ) # YES to exclude xchr data when calculating PCA
bld=$( awk '$1=="build" {print $2}' $paramfile )			# build of the dataset (18-19)
alpha=$( awk '$1=="alpha" {print $2}' $paramfile )			# alpha threshold to set
binaryplink=$( awk '$1=="plinkformat" {print $2}' $paramfile ) # ped = ped/map format, bed = binary plink format
minmaf=$( awk '$1=="maf" {print $2}' $paramfile )			# remove variants with MAF < minmaf
mindthresh=$( awk '$1=="missindiv" {print $2}' $paramfile )	# filter for missingness per individual
genothresh=$( awk '$1=="missgeno" {print $2}' $paramfile )	# filter for missingness per genotype
numpc=$( awk '$1=="numpc" {print $2}' $paramfile )			# number of PC's
pithresh=$( awk '$1=="related" {print $2}' $paramfile )		# threshold at which to remove related individuals
quant=$( awk '$1=="quant" {print $2}' $paramfile )			#is this a quantitative trait (either 0 or 1)
gender_check=$( awk '$1=="Gender_check" {print tolower($2)}' $paramfile )	# If set to NO, gender check will not be done
perform_pca=$( awk '$1=="perform_PCA" {print tolower($2)}' $paramfile )
HW_val=$( awk '$1=="HWE" {printf("%.20f\n", $2)}' $paramfile )	# The input needs to be scientific notation.
pval=$( awk '$1=="pval_cutoff" {print $2}' $paramfile )			# This value is used for X specific QC

if [ $eigstratloc != "." ]; then
	cp $eigstratloc/smartpca $eigstratloc/convertf ./
fi

if [ $exclind -eq 0 ]; then
	echo -n "" > ${fname}_exclind.remove
	cp ${fname}_exclind.remove ${fname}_male_exclind.remove
	cp ${fname}_exclind.remove ${fname}_female_exclind.remove
fi

totsnps=$(wc -l ${fname}.bim | awk '{print $1}')
bonf=$(echo "scale=20;$alpha/$totsnps" | bc)
bonf=$HW_val    # modified to take input from user
echo -e "HWE val is $bonf"

# Set PAR locations based on hg build
if [ $bld -eq 19 ]; then
	echo  '23 60001 2699520 par1\n' > pars.txt
	echo  '23 154931044 155260560 par2' >> pars.txt
fi

if [ $bld -eq 18 ]; then
	echo  '23 1 2709520 par1\n' > pars.txt
	echo  '23 154584238 154913754 par2' >> pars.txt
fi

if [ $bld -eq 17 ]; then
	echo  '23 1 2692881 par1\n' > pars.txt
	echo  '23 154494747 154824264 par2' >> pars.txt
fi

if [ $bld -eq 38 ]; then
	echo  '23 10001 2781479 par1\n' > pars.txt
	echo  '23 155701383 156030895 par2' >> pars.txt
fi

# remove pseudo-autosomal regions
echo -e "4: Removing pseudo-autosomal regions"
${xwasloc}/$xwasQ --bfile ${fname} --make-bed --out ${fname}_qc1 --exclude pars.txt --range

# check for wrong sex individuals (do automatically if > 1000 snps on the X chromosome) Otherwise do this manually
echo -e "5: Checking wrong-sex individuals"

xtot=`awk '$1==23 {print $0}' ./${fname}_qc1.bim | wc -l`
if [ $xtot -gt 999 ] && [ $gender_check = "yes" ]; then
	echo -e "performing the sex check"
	${xwasloc}/$xwasQ --bfile ${fname}_qc1 --check-sex --out ${fname}_psar_sexcheck
	cat ${fname}_psar_sexcheck.sexcheck | grep PROBLEM | awk '{print $1,$2,"sex"}' > ${fname}_sex_exclind.remove
else
	echo -e "Not performing the sex-check"
	touch ${fname}_psar_sexcheck
	touch ${fname}_sex_exclind.remove
fi

${xwasloc}/$xwasQ --bfile ${fname}_qc1 --make-bed --out ${fname}_qc2 --remove ${fname}_sex_exclind.remove
malenum=`awk   '$5 == 1' ${fname}_qc2.fam | wc -l`;
femalenum=`awk '$5 == 2' ${fname}_qc2.fam | wc -l`;

# Checking relatedness (IBD)
echo  '6 20000000 80000000 HLA' > hla.excl # generous buffer region for the HLA region
${xwasloc}/$xwasQ --bfile ${fname}_qc2 --indep-pairwise 50 5 0.3 --out ${fname}_qc2 --exclude hla.excl --range
allsnps=$(wc -l ${fname}_qc2.bim | awk '{print $1}')
thinamt=$(echo "scale=3;150000/$allsnps" | bc)

if [ $skipibd -eq 0 ]; then
	echo -e "6: Checking relatedness"
		
	## if the number of SNPs is less than 150,000, use all SNPs for PCA
	if [ $(echo $thinamt'>'1 | bc -l) -eq 1 ]; then
		echo "All SNPs used for IBD-based relatedness analysis";
		${xwasloc}/$xwasQ --bfile ${fname}_qc2 --genome --min ${pithresh} --extract ${fname}_qc2.prune.in --out ${fname}_rel
	else
		echo "150000 SNPs randomly chosen for IBD-based relatedness analysis"
		${xwasloc}/$xwasQ --bfile ${fname}_qc2 --thin $thinamt --genome --min ${pithresh} --extract ${fname}_qc2.prune.in --out ${fname}_rel
	fi
	
	perl ${xwasloc}/maximum_set_of_unrelated_individuals.pl   ${fname}_rel.genome  ${fname}_exclind.remove
	${xwasloc}/$xwasQ --bfile ${fname}_qc2 --make-bed --out ${fname}_qc3 --remove ${fname}_exclind.remove
else
	echo -e "6: Skipping IBD analysis"
	${xwasloc}/$xwasQ --bfile ${fname}_qc2 --make-bed --out ${fname}_qc3
fi
	
# PC
echo -e "7: Population stratification"

${xwasloc}/$xwasQ --bfile ${fname}_qc3 --extract ${fname}_qc2.prune.in --recode12 --out ${fname}_qc3
echo "RAn first xwas fine"
if [ $perform_pca = "yes" ]; then
	# convert to eigenstrat format
	${eigstratloc}/convertf.perl ${fname}_qc3
	# now get the PC's using eigenstrat
	${eigstratloc}/smartpca_loadings_mod.perl -i ${fname}_qc3.eigenstratgeno -a ${fname}_qc3.snp -b ${fname}_qc3.ind -k ${numpc} -o ${fname}_qc3.pca -p ${fname}_qc3.plot -e ${fname}_qc3.eval -l ${fname}_qc3.log -d ${fname}_qc3.load -m 5 -t 10 -s 6.0 -g ${xchrpca} > ${fname}_qc3.eig.log
	
	# Added by Lauren, catch if smartpca fails
	rc=$?
	if [ $rc -ne 0 ]; then
		echo "ERROR: smartpca exited with value $rc"
		echo "You may need to download the source code for smartpca and compile it for your machine."
		echo "You can download the source code from here: http://www.hsph.harvard.edu/alkes-price/software/"
		exit $rc
	fi

	# make the covariate file
	echo "Writing covariates to ${fname}.covar"
	awk 'BEGIN {OFS=" "} FNR==NR {f1[$2]=$1;next} ($1 in f1)  {$NF=""} {print f1[$1],$0}' ${fname}_qc2.fam ${fname}_qc3.pca.evec | grep -v eigvals: > ${fname}.covar
	# copying outlier individuals to the individuals to remove file
	cat ${fname}_qc3.log | grep REMOVED | awk '{print $3}' > temp.a
	echo -e "Writing the individuals to be excluded to ${fname}_exclind.remove"
	awk 'FNR==NR {f1[$0]; next} ($2 in f1) {print $1,$2,"pca"}' ./temp.a ${fname}_qc3.fam >> ${fname}_exclind.remove
fi

${xwasloc}/$xwasQ --bfile ${fname}_qc2 --remove ${fname}_exclind.remove --out ${fname}_qc4 --make-bed

# QC separately for male and female
if [ $malenum -eq 0 ] || [ $femalenum -eq 0 ]; then
	echo "Only one gender exist. Male # is $malenum; Female # is $femalenum"
	echo "QC is only performed in one gender"
	if [ $quant -eq 1 ]; then
		echo -e "8: Quality Control for Quantitative Traits"

		echo "8.1: HWE"
		${xwasloc}/$xwasQ --bfile ${fname}_qc4 --hardy --out ${fname}_1gender_hwe
		cat ${fname}_1gender_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_1gender_snp.exclude

		echo "8.2: MAF, missingness per SNP, missingness per individual"
		${xwasloc}/$xwasQ --bfile ${fname}_qc4 --make-bed --out ${fname}_qc5 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_1gender_snp.exclude
	else
		echo -e "8: Quality Control for Binary Traits"

		echo "8.1: HWE"
		${xwasloc}/$xwasQ --bfile ${fname}_qc4 --hardy --out ${fname}_1gender_hwe --filter-controls
		cat ${fname}_1gender_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_1gender_snp.exclude

		echo "8.2: Correlation between missingness and phenotype is omitted"
		${xwasloc}/$xwasQ --bfile ${fname}_qc4 --test-missing --out ${fname}_1gender_mcc
		# awk -v bf=$bonf '$5<bf {print $2}' ${fname}_1gender_mcc.missing >> ${fname}_1gender_snp.exclude

		echo "8.3: MAF, missingness per SNP, missingness per individual"
		${xwasloc}/$xwasQ --bfile ${fname}_qc4 --make-bed --out ${fname}_qc5 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_1gender_snp.exclude
	fi
else
	echo "QC is performed separately for each gender"
	if [ $quant -eq 1 ]; then
		echo -e "8: Quality Control for Quantitative Traits"

		echo "8.1: Quality Control for Male"
		${xwasloc}/$xwasQ --bfile ${fname}_qc4 --filter-males --make-bed --out ${fname}_male

		echo "8.1.1: HWE"
		${xwasloc}/$xwasQ --bfile ${fname}_male --hardy --out ${fname}_male_hwe
		cat ${fname}_male_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_male_snp.exclude

		echo "8.1.2: MAF, missingness per SNP, missingness per individual"
		${xwasloc}/$xwasQ --bfile ${fname}_male --make-bed --out ${fname}_male_qc5 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_male_snp.exclude

		echo "8.2: Quality Control for Female"
		${xwasloc}/$xwasQ --bfile ${fname}_qc4 --filter-females --make-bed --out ${fname}_female

		echo "8.2.1: HWE"
		${xwasloc}/$xwasQ --bfile ${fname}_female --hardy --out ${fname}_female_hwe
		cat ${fname}_female_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_female_snp.exclude

		echo "8.2.2: MAF, missingness per SNP, missingness per individual"
		${xwasloc}/$xwasQ --bfile ${fname}_female --make-bed --out ${fname}_female_qc5 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_female_snp.exclude

		echo -e "9: Merging male and female QC"
		perl ${xwasloc}/compare_SNPs_in_male_and_female.pl  ${fname}_female_qc5.bim  ${fname}_male_qc5.bim  ${fname}_sex_diff_snps.diff
		${xwasloc}/$xwasQ --bfile ${fname}_female_qc5 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_female_qc5_tmp
		${xwasloc}/$xwasQ --bfile ${fname}_male_qc5 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_male_qc5_tmp
		${xwasloc}/$xwasQ --bfile ${fname}_female_qc5_tmp --bmerge ${fname}_male_qc5_tmp.bed ${fname}_male_qc5_tmp.bim ${fname}_male_qc5_tmp.fam --make-bed --out ${fname}_final
	else
		echo -e "8: Quality Control for Binary Traits"

		echo "8.1: Quality Control for Male"
		${xwasloc}/$xwasQ --bfile ${fname}_qc4 --filter-males --make-bed --out ${fname}_male

		echo "8.1.1: HWE"
		${xwasloc}/$xwasQ --bfile ${fname}_male --hardy --out ${fname}_male_hwe --filter-controls
		cat ${fname}_male_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_male_snp.exclude

		echo "8.1.2: Correlation between missingness and phenotype is omitted"
		${xwasloc}/$xwasQ --bfile ${fname}_male --test-missing --out ${fname}_male_mcc
		# awk -v bf=$bonf '$5<bf {print $2}' ${fname}_male_mcc.missing >> ${fname}_male_snp.exclude

		echo "8.1.3: MAF, missingness per SNP, missingness per individual"
		${xwasloc}/$xwasQ --bfile ${fname}_male --make-bed --out ${fname}_male_qc5 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_male_snp.exclude

		echo "8.2: Quality Control for Female"
		${xwasloc}/$xwasQ --bfile ${fname}_qc4 --filter-females --make-bed --out ${fname}_female

		echo "8.2.1: HWE"
		${xwasloc}/$xwasQ --bfile ${fname}_female --hardy --out ${fname}_female_hwe --filter-controls
		cat ${fname}_female_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_female_snp.exclude

		echo "8.2.2: Correlation between missingness and phenotype is omitted"
		${xwasloc}/$xwasQ --bfile ${fname}_female --test-missing --out ${fname}_female_mcc
		# awk -v bf=$bonf '$5<bf {print $2}' ${fname}_female_mcc.missing >> ${fname}_female_snp.exclude

		echo "8.2.3: MAF, missingness per SNP, missingness per individual"
		${xwasloc}/$xwasQ --bfile ${fname}_female --make-bed --out ${fname}_female_qc5 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_female_snp.exclude

		echo -e "9: Merging male and female QC"
		perl ${xwasloc}/compare_SNPs_in_male_and_female.pl  ${fname}_female_qc5.bim  ${fname}_male_qc5.bim  ${fname}_sex_diff_snps.diff
		${xwasloc}/$xwasQ --bfile ${fname}_female_qc5 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_female_qc5_tmp
		${xwasloc}/$xwasQ --bfile ${fname}_male_qc5 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_male_qc5_tmp
		${xwasloc}/$xwasQ --bfile ${fname}_female_qc5_tmp --bmerge ${fname}_male_qc5_tmp.bed ${fname}_male_qc5_tmp.bim ${fname}_male_qc5_tmp.fam --make-bed --out ${fname}_final
	fi
fi

# Significant diff in MAF between males and females only for qualitative traits and only in controls
if [ $quant -eq 1 ]; then
	echo -e "10: Extracting X chromosome"
	${xwasloc}/$xwasQ --bfile ${fname}_final  --make-bed --out ${fname}_final_x  --chr X
	${xwasloc}/$xwasQ --bfile ${fname}_final_x --missing --out ${fname}_final_x_missing
else
	echo -e "10: Significant difference in MAF between males and females (only for binary traits and only in controls)"
	echo "Extracting X chromosome"
	# totx=$(awk '$1=="X" || $1==23 {print $0}' ${fname}_final.bim | wc -l)
	# bonfx=$(echo "scale=20;$alpha/$totx" | bc)
	# Take a parameter named pval cutoff and use it instead of bonfx
	bonfx=$pval
	${xwasloc}/$xwasQ --bfile ${fname}_final  --xwas --make-bed --out ${fname}_final_x  --chr X --freqdiff-x ${bonfx}
	${xwasloc}/$xwasQ --bfile ${fname}_final_x --missing --out ${fname}_final_x_missing
	# Run the --missing command and get the lmiss and imiss file
fi

if [ $eigstratloc != "." ]; then
	rm smartpca convertf
fi

# 5. Clean up after QC
if [ "$savemethod" = "none" ]; then # remove all intermediate files
	rm ${old_fname}.PC* ${fname}.bed ${fname}.bim ${fname}.fam
	rm ${fname}_exclind.remove ${fname}_female* ${fname}_male*
	rm ${fname}_psar* ${fname}_qc* 
	rm ${fname}_sex* ${old_fname}_snp_only.Sib* pars.txt temp.a hla.excl temp_params.txt
	if [ $skipibd -eq 0 ]; then
		rm ${fname}_rel*
	fi
	rm ${old_fname}_snp_only*
elif [ "$savemethod" = "logs" ]; then # save logs, remove other intermediates
	LOGS=${old_fname}_QC_logs
	mkdir -p ${LOGS}
	mv temp_params.txt *.log $LOGS/
	rm ${old_fname}.PC* ${fname}.bed ${fname}.bim ${fname}.fam
	rm ${fname}_exclind.remove ${fname}_female* ${fname}_male*
	rm ${fname}_psar* ${fname}_qc*
	rm ${fname}_sex* ${old_fname}_snp_only.Sib* pars.txt temp.a hla.excl
	if [ $skipibd -eq 0 ]; then
		rm ${fname}_rel*
	fi
	rm ${old_fname}_snp_only*
else # keep everything
	INTER=${old_fname}_QC_intermediate_files
	mkdir -p ${INTER}
	mv ${old_fname}.PC* ${fname}.bed ${fname}.bim ${fname}.fam $INTER/
	mv ${fname}_exclind.remove ${fname}_female* ${fname}_male* $INTER/
	mv ${fname}_psar* ${fname}_qc* $INTER/
	mv ${fname}_sex* ${old_fname}_snp_only.Sib* pars.txt temp.a hla.excl temp_params.txt $INTER/
	if [ $skipibd -eq 0 ]; then
		mv ${fname}_rel* $INTER/
	fi
	mv ${old_fname}_snp_only* $INTER/
fi

echo -e "Performing the frequency calculation using Plink within xwas executable\n"
parameter_file="$key"
echo -e "Reading params from ${parameter_file}\n"
af_input_file=${fname}_final_x	# input file for allele freq test, Performing analysis on X Chr output
default_case_control=$( awk '$1=="default_case_control" {print tolower($2)}' ${parameter_file} )
cases=$( awk '$1=="cases" {print tolower($2)}' ${parameter_file} ) # YES to perform them separately
controls=$( awk '$1=="controls" {print tolower($2)}' ${parameter_file} )
default_male_female=$( awk '$1=="default_male_female" {print tolower($2)}' ${parameter_file})
males=$( awk '$1=="males" {print tolower($2)}' ${parameter_file} ) # YES to perform them separately
females=$( awk '$1=="females" {print tolower($2)}' ${parameter_file} )
p_comb_method=$( awk '$1=="P_comb_method" {print tolower($2)}' ${parameter_file} )
ci=$( awk '$1=="confidence_interval" {print $2}' ${parameter_file} )
outfile_name=$(echo "${fname}_final_x" | sed 's/\.preprocessed_final_x//' )

declare -A CASE_CONTROL=( [default]="" [cases]="--filter-cases" [controls]="--filter-controls" )
declare -A MALE_FEMALE=( [default]="" [male]="--filter-males" [female]="--filter-females")

if [ $default_case_control = "no" ]
then
	unset CASE_CONTROL[default]
fi

if [ $cases = "no" ]
then
	unset CASE_CONTROL[cases]
fi

if [ $controls = "no" ]
then
	unset CASE_CONTROL[controls]
fi

if [ $default_male_female = "no" ]
then
	unset MALE_FEMALE[default]
fi

if [ $males = "no" ]
then
	unset MALE_FEMALE[male]
fi

if [ $females = "no" ]
then
	unset MALE_FEMALE[female]
fi
 
for values in "${!CASE_CONTROL[@]}"
do
	for gender_values in "${!MALE_FEMALE[@]}"
	do
		${xwasloc}/xwas --bfile $af_input_file --noweb --freq --out ${outfile_name}_${gender_values}_${values} ${CASE_CONTROL[${values}]} ${MALE_FEMALE[${gender_values}]}
	done
done

# Sex stratified association test and sex-difference test
if [ $p_comb_method = "stouffers" ]
then
	method="--stouffers"
else
	method="--fishers"
fi

${xwasloc}/xwas --bfile $af_input_file --noweb --xwas --sex-diff --strat-sex $method --ci $ci --gc --out ${outfile_name}_sex_diff_strat

# Analysing the output of the sex-stratified test and allele frequency steps
combined_pvalue_column=$(awk '{for (i=1; i<=NF; i++) { if (tolower($i) ~ /p_comb_fisher|p_comb_stouffer/) print i}; exit}' ${outfile_name}_sex_diff_strat.xstrat.logistic | head -n1)
awk -v col=$combined_pvalue_column -v p=$pval '{if (NR>1 && $col < p) print $2}' ${outfile_name}_sex_diff_strat.xstrat.logistic > ${outfile_name}_SNP.significant

if [ -f "${outfile_name}_male_cases.frq" ] && [ -f "${outfile_name}_male_controls.frq" ]; then
	grep -w -f ${outfile_name}_SNP.significant ${outfile_name}_male_cases.frq | awk 'BEGIN{OFS="\t"} {print $2,$5}' > ${outfile_name}_mcase.dat
	grep -w -f ${outfile_name}_SNP.significant ${outfile_name}_male_controls.frq | awk '{print $5}' > ${outfile_name}_mcontrol.dat
	touch ${outfile_name}_male_case_control_SNP.significant
	echo -e "SNP\tmale_cases\tmale_controls" > ${outfile_name}_male_case_control_SNP.significant
	paste ${outfile_name}_mcase.dat ${outfile_name}_mcontrol.dat >> ${outfile_name}_male_case_control_SNP.significant
	rm ${outfile_name}_mcase.dat ${outfile_name}_mcontrol.dat
else
	echo -e "male case/ control files are missing"
fi

if [ -f "${outfile_name}_female_cases.frq" ] && [ -f "${outfile_name}_female_controls.frq" ]; then
	grep -w -f ${outfile_name}_SNP.significant ${outfile_name}_female_cases.frq | awk 'BEGIN{OFS="\t"} {print $2,$5}' > ${outfile_name}_fcase.dat
	grep -w -f ${outfile_name}_SNP.significant ${outfile_name}_female_controls.frq | awk '{print $5}' > ${outfile_name}_fcontrol.dat
	touch ${outfile_name}_female_case_control_SNP.significant
	echo -e "SNP\tfemale_cases\tfemale_controls" > ${outfile_name}_female_case_control_SNP.significant
	paste ${outfile_name}_fcase.dat ${outfile_name}_fcontrol.dat >> ${outfile_name}_female_case_control_SNP.significant
	rm ${outfile_name}_fcase.dat ${outfile_name}_fcontrol.dat
else
	echo -e "female case/control files are missing"
fi