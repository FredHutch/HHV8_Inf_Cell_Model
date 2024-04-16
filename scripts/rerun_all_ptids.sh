#!/bin/bash
rm -f single_runs.log
model=$1
shift
for i in "$@"
do
    ptid="ptid_""$i"
    echo "For $ptid ..." >> single_runs.log
    if [ ! -d $ptid\_$model ]
    then
	mkdir $ptid\_$model
    fi
    cd $ptid\_$model
    cp ../$ptid\_$model.in .
    cp ../Test.ttf .
    echo "N_runs 10" >> $ptid\_$model.in
    echo "beta_std 0" >> $ptid\_$model.in
    echo "fpos_std 0" >> $ptid\_$model.in
    echo "an_std 0" >> $ptid\_$model.in
    echo "betae_std 0" >> $ptid\_$model.in
    echo "log_p_std 0" >> $ptid\_$model.in
    echo "latent_inf_std 0" >> $ptid\_$model.in
    echo "r_std 0" >> $ptid\_$model.in
    echo "exp_days_std 0" >> $ptid\_$model.in
    echo "alpha_std 0" >> $ptid\_$model.in
    echo "kappa_std 0" >> $ptid\_$model.in

    ../run_pdf_set.pl 1 $ptid\_$model.in ../high_cohort.crit
    ../distrib_runs.pl hhv8_sim.dat1 > ../$ptid\_$model.dat1.csv
    cd ..
done

