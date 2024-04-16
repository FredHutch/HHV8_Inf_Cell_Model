#!/bin/sh
prefix=$1
ptid=$2
model=$3
txtfile=$4

./indiv_rollup.pl $txtfile  $prefix.crit $prefix\_summary.crit $ptid > $prefix\_$model.csv
head -n 1 $prefix\_$model.csv > top100_$prefix\_$model.csv

# sort by shed/med hhv8 score
grep -v Run $prefix\_$model.csv > temp_scores.csv
sort -n -k82  -t',' temp_scores.csv |head -n 100 >>  top100_$prefix\_$model.csv
./top100_dist.pl top100_$prefix\_$model.csv $model >> $prefix\_$model.in
