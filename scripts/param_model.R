#!/bin/bash
model=$1
shift
for i in "$@"
do
echo "For ptid $i, model $model..."
#Run,mean beta,mean log_p,mean latent_inf,mean theta,mean delta,mean r,mean c,mean rho,mean inf,mean rinf,mean eclipse,mean hill,mean To,mean fpos,mean an,regs,avg shedrate,score1,score2,score3,score,sim1_1,act1_1,sim1_2,act1_2,sim1_3,act1_3,sim1_4,act1_4,sim1_5,act1_5,new_score1,sim2_0,act2_0,sim2_1,act2_1,sim2_2,act2_2,sim2_3,act2_3,sim2_4,act2_4,new_score2,sim3_0,act3_0,sim3_1,act3_1,sim3_2,act3_2,sim3_3,act3_3,sim3_4,act3_4,sim3_5,act3_5,new_score3,new_score,AIC_k,AIC_score,mean kappa,mean cd8_ic50,total new score,old score,shedding,shedding target,shedding score,med log hhv8,med log hhv8 target,med hhv8 score,var log hhv8, var target,var hhv8 score,ptid,exp days,alpha,avg ibirths,avg tbirths,peak log hhv8, peak target,peak hhv8 score,three score,betae
cut -f2,3,4,15,83 -d',' top100_ptid_$i\_$model.csv | head -n 11|awk -F',' '{if (NR > 1){n++;beta+=$1;fpos+=$2;log_p+=$3;lat_react+=$4;betae+=$5;}}END{printf("%d params: beta= %e fpos= %g lat_react= %g log_p= %g betae=%e\n",n,beta/n, fpos/n,lat_react/n,log_p/n,betae/n);}'
echo
done
