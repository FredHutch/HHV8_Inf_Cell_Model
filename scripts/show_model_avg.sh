#!/bin/sh
model=$1
file=$2
echo -n "For model $model: "
cat $file | awk '{if ($2=="scores:") {printf("%g,%g,%g,%g\n",$4,$8,$10,$12);}}'|awk -F',' '{n++;sum+=$1;sum1+=$2;sum2+=$3;sum3+=$4;}END{printf("%d scores, %g avg, %g shed, %g hhv8, %g peak\n",n,sum/n,sum1/n,sum2/n,sum3/n);}'
