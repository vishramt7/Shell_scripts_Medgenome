#!/bin/bash
input="dat1.txt"
j=0 # introduce a variable
while read line; do
    ((j++))
    echo "${j} ${line}" > line_no_${j}.txt # Saving it to a new file
    wc -l line_no_${j}.txt > second_command_output_${j}.txt # Open the file saved above and use it as input
done < "$input"