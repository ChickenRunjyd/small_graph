#!/bin/bash
# File: diff.bash
# Date: 15:25 10/30/2015
# Description: diff the output of CPU and GPU BC

file_1="bc_bottom_up_gpu."
file_2="bc_bottom_up_cpu_debug."
deal () {
    diff $file_1$1 $file_2$1 > diff.$1
}

argu[0]="dist"
argu[1]="sa"
argu[2]="sp_count"
argu[3]="bc"

for index in 0 1 2 3
do
    deal ${argu[$index]}
done
