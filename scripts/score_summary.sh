#!/bin/bash
model=$1
shift
for i in "$@"
do
echo "For ptid $i..."
cut -f82 -d',' top100_ptid_$i\_$model.csv | head -n 11|awk 'BEGIN{mins=1000;maxs=0;}{if (NR > 1){n++;sum+=$1;if (mins > $1){mins=$1;}if (maxs < $1){maxs=$1;}}}END{printf("model 1, %d scores: Avg=%g (max=%g, min=%g)\n",n,sum/n, maxs,mins);}'
echo
done
