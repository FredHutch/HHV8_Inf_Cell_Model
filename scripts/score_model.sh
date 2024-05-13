#!/bin/bash
model=$1
shift
for i in "$@"
do
echo "For ptid $i, model $model..."
cut -f67,70,81,82 -d',' top100_ptid_$i\_$model.csv | head -n 11|awk -F',' 'BEGIN{mins=1000;maxs=0;}{if (NR > 1){n++;sum+=$4;sum1+=$1;sum2+=$2;sum3+=$3;if (mins > $4){mins=$4;}if (maxs < $4){maxs=$4;}}}END{printf("%d scores: Avg= %g (max=%g, min=%g) shed= %g hhv8= %g peak= %g\n",n,sum/n, maxs,mins,sum1/n,sum2/n,sum3/n);}'
echo
done
