#!/bin/bash
# File: debug.bash
# Date: 18:52 10/31/2015
# Description: diff the output of CPU and GPU BC
#               from every vertex to others

path_cpu="/home/yuede/small_graph/bc_debug/";
path_gpu="/home/yuede/small_graph/bc/"
path_result="/home/yuede/small_graph/result_bc/"
diff_bc="diff.bc"
diff_dist="diff.dist"
diff_sa="diff.sa"
diff_sp_count="diff.sp_count"

deal () {

    cd $path_cpu
    #start_vert=$1 end_vert=$(($1+1)) make test
    start_vert=$1 end_vert=$2 make test_debug
    
    cd $path_gpu
    start_vert=$1 end_vert=$2 make test_debug

    cd $path_result
    ./diff.bash

    #actualsize=$(wc -c <"$file")
    if [ -s diff.dist ]
    then
        echo $1 diff.dist wrong 
        exit 1
    fi

    if [ -s diff.sa ]
    then
        echo $1 diff.sa wrong
        #exit 1
    fi
    
    if [ -s diff.sp_count ]
    then
        echo $1 diff.sp_count wrong 
        #exit 1
    fi
    
    if [ -s diff.bc ]
    then
        echo $1 diff.bc wrong 
        exit 1
    fi
}
#cd $path_cpu

#deal 0 218

#for index in `seq $1 218`;
for index in `seq 0 $((218-$1))`;
do
    #echo $index $(($index+2))
    #echo $index $(($index+$1))
    #deal $index $(($index 
    deal $index $(($index+$1))
done
