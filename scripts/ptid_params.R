library(RColorBrewer)
args<-commandArgs(trailingOnly=T)
ptid<-args[1]

ptid_data = read.csv(paste0("top100_",ptid,".csv"))
#Run,mean beta,mean log_p,mean latent_inf,mean theta,mean delta,mean r,mean c,mean rho,mean inf,mean rinf,mean eclipse,mean betae,beta * p,mean To,mean fpos,mean an,regs,avg shedrate,score1,score2,score3,score,sim1_1,act1_1,sim1_2,act1_2,sim1_3,act1_3,sim1_4,act1_4,sim1_5,act1_5,new_score1,sim2_0,act2_0,sim2_1,act2_1,sim2_2,act2_2,sim2_3,act2_3,sim2_4,act2_4,new_score2,sim3_0,act3_0,sim3_1,act3_1,sim3_2,act3_2,sim3_3,act3_3,sim3_4,act3_4,sim3_5,act3_5,new_score3,new_score,AIC_k,AIC_score,mean kappa,mean cd8_ic50,total new score,shedding,shedding target,med log hhv8,med log hhv8 target,peak log hhv8, peak target,peak hhv8 score,ptid,mean To

ptid_data$Run
ptid_data$mean.beta
ptid_data$mean.log_p
ptid_data$mean.latent_inf

min(ptid_data$mean.latent_inf)
max(ptid_data$mean.latent_inf)
min_li=log10(min(ptid_data$mean.latent_inf))
max_li=log10(max(ptid_data$mean.latent_inf))
min_fp=log10(min(ptid_data$mean.fpos))
max_fp=log10(max(ptid_data$mean.fpos))
min_lb=log10(min(ptid_data$mean.beta))
max_lb=log10(max(ptid_data$mean.beta))

low_log_beta=-11
high_log_beta=-9

low_log_beta_p=low_log_beta+4
high_log_beta_p=high_log_beta+4

low_log_fpos=-3
high_log_fpos=0

low_log_latent_inf=-2
high_log_latent_inf=0

pdf(paste0(ptid,"_misc_params.pdf"),height=8,width=8)
par(mfrow = c(2,2), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

# latent_reacts
y_lim=c(min_li,max_li)
print(y_lim)

latent_infs=log10(ptid_data$mean.latent_inf)
boxplot(latent_infs, range=0.0, horizontal=FALSE, col=c("red", "green"),
	main="Log Activation rate",
        ylim=y_lim,xlab = "", ylab = "log latent inf")


# latent_inf
x_lim=c(min(low_log_latent_inf,min_li),
	max(high_log_latent_inf,max_li))
y_lim=c(0,max(ptid_data$total.new.score))

plot(x = x_lim, y = y_lim, type = "n", las = 1, xlab = "latent inf", ylab = "score",
       xlim = x_lim, ylim = y_lim)
points(log10(ptid_data$mean.latent_inf), ptid_data$total.new.score, col = "red", pch = 1)

# fpos
y_lim=c(min_fp,max_fp)
print(y_lim)

fpos=log10(as.numeric(ptid_data$mean.fpos))
boxplot(fpos, range=0.0, horizontal=FALSE, col=c("red", "green"),
	main="Log Killing Rate",
        ylim=y_lim,xlab = "", ylab = "log fpos")


# fpos
x_lim=c(min(low_log_fpos,log10(as.numeric(ptid_data$mean.fpos))),
	max(high_log_fpos,log10(as.numeric(ptid_data$mean.fpos))))
y_lim=c(0,max(ptid_data$total.new.score))

plot(x = x_lim, y = y_lim, type = "n", las = 1, xlab = "log fpos", ylab = "score",
       xlim = x_lim, ylim = y_lim)
points(log10(ptid_data$mean.fpos), ptid_data$total.new.score, col = "red", pch = 1)

dev.off()

pdf(paste0(ptid,"_viral_params.pdf"),height=8,width=8)
par(mfrow = c(2,2), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

# beta
y_lim=c(min_lb,max_lb)
print(y_lim)

beta=log10(as.numeric(ptid_data$mean.beta))
boxplot(beta, range=0.0, horizontal=FALSE, col=c("red", "green"),
	main="Log Infectivity",
        ylim=y_lim,xlab = "", ylab = "log beta")


# beta
x_lim=c(min(low_log_beta,log10(as.numeric(ptid_data$mean.beta))),
	max(high_log_beta,log10(as.numeric(ptid_data$mean.beta))))
y_lim=c(0,max(ptid_data$total.new.score))

plot(x = x_lim, y = y_lim, type = "n", las = 1, xlab = "log beta", ylab = "score",
       xlim = x_lim, ylim = y_lim)
points(log10(ptid_data$mean.beta), ptid_data$total.new.score, col = "red", pch = 1)

dev.off()

