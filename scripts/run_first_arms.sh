#!/bin/bash
if [ "$1" != "" ]
then
model=$1
else
model=1
fi

if [ "$2" != "" ]
then
sets=$2
else
sets=10
fi

if [ "$3" != "" ]
then
starts=$3
else
starts=1
fi

if [ "$4" != "" ]
then
count2=$4
else
count2=1
fi

critfile="../no_betae_high.crit"

infile="arms_"
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
    cp ../Test.ttf .
    cp ../$infile .
    rm -f slurm*
    sbatch -n 1 -t 20:00:00 ../run_first_sets.pl $count2 $infile $critfile $model
    cd ..

    sets=`expr $sets - 1`
    count=`expr $count + 1`
done
