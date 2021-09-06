#!/bin/bash

declare -A assArray2=( [HDD]=Samsung [Monitor]=Dell [Keyboard]=A4Tech )
echo ${assArray2[HDD]}
declare -A CASE_CONTROL=( [default]="" [cases]="--filter-cases" [controls]="--filter-controls" )
declare -A MALE_FEMALE=( [default]="" [male]="--filter-males" [female]="--filter-females")

unset CASE_CONTROL[cases]
for keys in "${!CASE_CONTROL[@]}"
do
    for gender_keys in "${!MALE_FEMALE[@]}"
    do
        echo "${keys} => ${CASE_CONTROL[$keys]} , ${gender_keys} => ${MALE_FEMALE[$gender_keys]}"
    done
done
