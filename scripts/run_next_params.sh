#!/bin/bash
prefix=$1
if [ "$2" != "" ]
then
model=$2
else
model=1
fi

if [ "$3" != "" ]
then
sets=$3
else
sets=10
fi

if [ "$4" != "" ]
then
starts=$4
else
starts=1
fi

if [ "$5" != "" ]
then
count2=$5
else
count2=1
fi

critfile="../../data/high_shedders.crit"

infile=$prefix
infile+="_"
infile+=$model
infile+=".in"

# run all sims varying these parameters (via run_param_sets.pl)

count=$starts
while [ $sets -gt 0 ]
do
    thisdir="set_"
    thisdir+=$count
    if [ ! -d $thisdir ]
    then
    mkdir $thisdir
    fi
    cd $thisdir
    cp ../$prefix\_$model.in $prefix\_$model.in
    rm -f slurm*
    sbatch -n 1 -t 20:00:00 ../../scripts/run_single_set.pl $count2 $infile $critfile $model
    cd ..

    sets=`expr $sets - 1`
    count=`expr $count + 1`
done
