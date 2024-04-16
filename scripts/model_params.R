library(RColorBrewer)
args<-commandArgs(trailingOnly=T)
model<-args[1]
tag<-args[2]

ptid_data = read.csv(paste0("params_",model,"_",tag,".csv"))
#ptid,an,beta,log_p,latent_inf,r,exp_days,betae

low_beta=1
high_beta=20

low_an=0;
high_an=3;

low_log_p=3
high_log_p=6

low_log_fpos=-2
high_log_fpos=2

low_log_latent_inf=-2
high_log_latent_inf=1

low_exp_days=7;
high_exp_days=20;

low_log_r=1;
high_log_r=3;

if (model==5) {
    pdf(paste0("model_",model,"_",tag,"_misc_params.pdf"),height=6,width=12)
    par(mfrow = c(1,3), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))
} else {
    pdf(paste0("model_",model,"_",tag,"_misc_params.pdf"),height=6,width=8)
    par(mfrow = c(1,2), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))
}

# latent_reacts
y_lim=c(low_log_latent_inf,high_log_latent_inf)
print(y_lim)

latent_infs=log10(ptid_data$latent_inf)
boxplot(latent_infs, range=0.0, horizontal=FALSE, col=c("red"),
	main="Log Activation rate",
        ylim=y_lim,xlab = "", ylab = "log latent inf")

# death rate
y_lim=c(low_an,high_an)
print(y_lim)

boxplot(ptid_data$an, range=0.0, horizontal=FALSE, col=c("green"),
	main="Inf Cell Death rate",
        ylim=y_lim,xlab = "", ylab = "death rate")

if (model==5) {
    # exp_days
    y_lim=c(low_exp_days,high_exp_days)
    print(y_lim)

    boxplot(ptid_data$exp_days, range=0.0, horizontal=FALSE, col=c("green"),
	    main="Expansion Days",
	    ylim=y_lim,xlab = "", ylab = "exp_days")
}

dev.off()

pdf(paste0("model_",model,"_",tag,"_tcell_params.pdf"),height=6,width=8)
par(mfrow = c(1,2), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

# fpos
y_lim=c(low_log_fpos,high_log_fpos)
print(y_lim)

fpos=log10(as.numeric(ptid_data$fpos))
boxplot(fpos, range=0.0, horizontal=FALSE, col=c("red"),
	main="Log Killing Rate",
        ylim=y_lim,xlab = "", ylab = "log fpos")

# r
y_lim=c(low_log_r,high_log_r)
print(y_lim)

r=log10(as.numeric(ptid_data$r))
boxplot(r, range=0.0, horizontal=FALSE, col=c("green"),
	main="Log Inf Cell IC50",
        ylim=y_lim,xlab = "", ylab = "log r")

dev.off()

pdf(paste0("model_",model,"_",tag,"_viral_params.pdf"),height=6,width=8)
par(mfrow = c(1,2), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

# beta
y_lim=c(low_beta, high_beta)
print(y_lim)

boxplot(ptid_data$beta, range=0.0, horizontal=FALSE, col=c("red"),
	main="Infectivity",
        ylim=y_lim,xlab = "", ylab = "beta")

# log_p
y_lim=c(low_log_p,high_log_p)
print(y_lim)

boxplot(ptid_data$log_p, range=0.0, horizontal=FALSE, col=c("green"),
	main="Log Viral Production",
        ylim=y_lim,xlab = "", ylab = "log p")


dev.off()

