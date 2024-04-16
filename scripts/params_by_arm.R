library(RColorBrewer)
# Load ggplot2
library(ggplot2)
library(tidyverse)

args<-commandArgs(trailingOnly=T)
model<-args[1]

both_pos = read.csv(paste0("both_pos_",model,"_sim.csv"))
both_neg = read.csv(paste0("both_neg_",model,"_sim.csv"))
hiv_pos = read.csv(paste0("hiv_pos_",model,"_sim.csv"))
ks_pos = read.csv(paste0("ks_pos_",model,"_sim.csv"))
#ptid,an,beta,log_p,latent_inf,r,exp_days,betae

shedding=c(both_pos$shedrate,both_neg$shedrate,hiv_pos$shedrate,ks_pos$shedrate)
log_hhv8=c(both_pos$median.log.hhv8,both_neg$median.log.hhv8,hiv_pos$median.log.hhv8,ks_pos$median.log.hhv8)
arms=c(rep("HIV+/KS+",nrow(both_pos)),rep("HIV-/KS-",nrow(both_neg)),rep("HIV+/KS-",nrow(hiv_pos)),rep("HIV-/KS+",nrow(ks_pos)))
arm_colors=c(rep("red",nrow(both_pos)),rep("green",nrow(both_neg)),rep("orange",nrow(hiv_pos)),rep("blue",nrow(ks_pos)))

#low_beta=1
#high_beta=20
low_beta=min(both_pos$beta,both_neg$beta,hiv_pos$beta,ks_pos$beta)
high_beta=max(both_pos$beta,both_neg$beta,hiv_pos$beta,ks_pos$beta)

mean1=mean(both_pos$beta)
sd1=sd(both_pos$beta)

mean2=mean(both_neg$beta)
sd2=sd(both_neg$beta)

mean3=mean(hiv_pos$beta)
sd3=sd(hiv_pos$beta)

mean4=mean(ks_pos$beta)
sd4=sd(ks_pos$beta)
print("beta:")
print(paste0("both_pos (mean=",formatC(mean1,digits=3,format="f"),",sd=",formatC(sd1,digits=3,format="f"),")"))
print(paste0("both_neg (mean=",formatC(mean2,digits=3,format="f"),",sd=",formatC(sd2,digits=3,format="f"),")"))
print(paste0("hiv_pos (mean=",formatC(mean3,digits=3,format="f"),",sd=",formatC(sd3,digits=3,format="f"),")"))
print(paste0("ks_pos (mean=",formatC(mean4,digits=3,format="f"),",sd=",formatC(sd4,digits=3,format="f"),")"))

#low_an=0;
#high_an=3;
low_an=min(both_pos$an,both_neg$an,hiv_pos$an,ks_pos$an)
high_an=max(both_pos$an,both_neg$an,hiv_pos$an,ks_pos$an)

mean1=mean(both_pos$an)
sd1=sd(both_pos$an)

mean2=mean(both_neg$an)
sd2=sd(both_neg$an)

mean3=mean(hiv_pos$an)
sd3=sd(hiv_pos$an)

mean4=mean(ks_pos$an)
sd4=sd(ks_pos$an)
print(paste0("an:"))
print(paste0("both_pos (mean=",formatC(mean1,digits=3,format="f"),",sd=",formatC(sd1,digits=3,format="f"),")"))
print(paste0("both_neg (mean=",formatC(mean2,digits=3,format="f"),",sd=",formatC(sd2,digits=3,format="f"),")"))
print(paste0("hiv_pos (mean=",formatC(mean3,digits=3,format="f"),",sd=",formatC(sd3,digits=3,format="f"),")"))
print(paste0("ks_pos (mean=",formatC(mean4,digits=3,format="f"),",sd=",formatC(sd4,digits=3,format="f"),")"))

#low_log_p=3
#high_log_p=8
low_log_p=min(both_pos$log_p,both_neg$log_p,hiv_pos$log_p,ks_pos$log_p)
high_log_p=max(both_pos$log_p,both_neg$log_p,hiv_pos$log_p,ks_pos$log_p)

mean1=mean(both_pos$log_p)
sd1=sd(both_pos$log_p)

mean2=mean(both_neg$log_p)
sd2=sd(both_neg$log_p)

mean3=mean(hiv_pos$log_p)
sd3=sd(hiv_pos$log_p)

mean4=mean(ks_pos$log_p)
sd4=sd(ks_pos$log_p)
print(paste0("log_p:"))
print(paste0("both_pos (mean=",formatC(mean1,digits=3,format="f"),",sd=",formatC(sd1,digits=3,format="f"),")"))
print(paste0("both_neg (mean=",formatC(mean2,digits=3,format="f"),",sd=",formatC(sd2,digits=3,format="f"),")"))
print(paste0("hiv_pos (mean=",formatC(mean3,digits=3,format="f"),",sd=",formatC(sd3,digits=3,format="f"),")"))
print(paste0("ks_pos (mean=",formatC(mean4,digits=3,format="f"),",sd=",formatC(sd4,digits=3,format="f"),")"))

#low_log_fpos=-2
#high_log_fpos=2
low_log_fpos=min(log10(both_pos$fpos),log10(both_neg$fpos),log10(hiv_pos$fpos),log10(ks_pos$fpos))
high_log_fpos=max(log10(both_pos$fpos),log10(both_neg$fpos),log10(hiv_pos$fpos),log10(ks_pos$fpos))

mean1=mean(both_pos$fpos)
sd1=sd(both_pos$fpos)

mean2=mean(both_neg$fpos)
sd2=sd(both_neg$fpos)

mean3=mean(hiv_pos$fpos)
sd3=sd(hiv_pos$fpos)

mean4=mean(ks_pos$fpos)
sd4=sd(ks_pos$fpos)
print(paste0("fpos:"))
print(paste0("both_pos (mean=",formatC(mean1,digits=3,format="f"),",sd=",formatC(sd1,digits=3,format="f"),")"))
print(paste0("both_neg (mean=",formatC(mean2,digits=3,format="f"),",sd=",formatC(sd2,digits=3,format="f"),")"))
print(paste0("hiv_pos (mean=",formatC(mean3,digits=3,format="f"),",sd=",formatC(sd3,digits=3,format="f"),")"))
print(paste0("ks_pos (mean=",formatC(mean4,digits=3,format="f"),",sd=",formatC(sd4,digits=3,format="f"),")"))

#low_log_latent_inf=-2
#high_log_latent_inf=2
low_log_latent_inf=min(log10(both_pos$latent_inf),log10(both_neg$latent_inf),log10(hiv_pos$latent_inf),log10(ks_pos$latent_inf))
high_log_latent_inf=max(log10(both_pos$latent_inf),log10(both_neg$latent_inf),log10(hiv_pos$latent_inf),log10(ks_pos$latent_inf))

#low_alpha=0.1;
#high_alpha=0.5;
low_alpha=min(both_pos$alpha,both_neg$alpha,hiv_pos$alpha,ks_pos$alpha)
high_alpha=max(both_pos$alpha,both_neg$alpha,hiv_pos$alpha,ks_pos$alpha)

mean1=mean(both_pos$alpha)
sd1=sd(both_pos$alpha)

mean2=mean(both_neg$alpha)
sd2=sd(both_neg$alpha)

mean3=mean(hiv_pos$alpha)
sd3=sd(hiv_pos$alpha)

mean4=mean(ks_pos$alpha)
sd4=sd(ks_pos$alpha)
print(paste0("alpha:"))
print(paste0("both_pos (mean=",formatC(mean1,digits=3,format="f"),",sd=",formatC(sd1,digits=3,format="f"),")"))
print(paste0("both_neg (mean=",formatC(mean2,digits=3,format="f"),",sd=",formatC(sd2,digits=3,format="f"),")"))
print(paste0("hiv_pos (mean=",formatC(mean3,digits=3,format="f"),",sd=",formatC(sd3,digits=3,format="f"),")"))
print(paste0("ks_pos (mean=",formatC(mean4,digits=3,format="f"),",sd=",formatC(sd4,digits=3,format="f"),")"))

#low_log_r=2;
#high_log_r=3;
low_log_r=min(log10(both_pos$r),log10(both_neg$r),log10(hiv_pos$r),log10(ks_pos$r))
high_log_r=max(log10(both_pos$r),log10(both_neg$r),log10(hiv_pos$r),log10(ks_pos$r))

mean1=mean(both_pos$r)
sd1=sd(both_pos$r)

mean2=mean(both_neg$r)
sd2=sd(both_neg$r)

mean3=mean(hiv_pos$r)
sd3=sd(hiv_pos$r)

mean4=mean(ks_pos$r)
sd4=sd(ks_pos$r)
print(paste0("r:"))
print(paste0("both_pos (mean=",formatC(mean1,digits=3,format="f"),",sd=",formatC(sd1,digits=3,format="f"),")"))
print(paste0("both_neg (mean=",formatC(mean2,digits=3,format="f"),",sd=",formatC(sd2,digits=3,format="f"),")"))
print(paste0("hiv_pos (mean=",formatC(mean3,digits=3,format="f"),",sd=",formatC(sd3,digits=3,format="f"),")"))
print(paste0("ks_pos (mean=",formatC(mean4,digits=3,format="f"),",sd=",formatC(sd4,digits=3,format="f"),")"))

# latent_reacts
y_lim=c(low_log_latent_inf,high_log_latent_inf)
latent_infs=list(log10(both_pos$latent_inf),log10(both_neg$latent_inf),log10(hiv_pos$latent_inf),log10(ks_pos$latent_inf))

mean1=mean(both_pos$latent_inf)
sd1=sd(both_pos$latent_inf)

mean2=mean(both_neg$latent_inf)
sd2=sd(both_neg$latent_inf)

mean3=mean(hiv_pos$latent_inf)
sd3=sd(hiv_pos$latent_inf)

mean4=mean(ks_pos$latent_inf)
sd4=sd(ks_pos$latent_inf)
print(paste0("latent_inf:"))
print(paste0("both_pos (mean=",formatC(mean1,digits=3,format="f"),",sd=",formatC(sd1,digits=3,format="f"),")"))
print(paste0("both_neg (mean=",formatC(mean2,digits=3,format="f"),",sd=",formatC(sd2,digits=3,format="f"),")"))
print(paste0("hiv_pos (mean=",formatC(mean3,digits=3,format="f"),",sd=",formatC(sd3,digits=3,format="f"),")"))
print(paste0("ks_pos (mean=",formatC(mean4,digits=3,format="f"),",sd=",formatC(sd4,digits=3,format="f"),")"))

pdf(paste0("arms_",model,"_avgs_latent_inf.pdf"),height=6,width=6)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

boxplot(latent_infs, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Log reactivation rate (all participants)" ,
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
        ylim=y_lim,xlab = "", ylab = "log reactivation rate")
#axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

legend("bottomleft", inset=.02, 
       c("HIV+/KS+", "HIV-/KS-", "HIV+/KS-", "HIV-/KS+"), 
       fill= c("red", "green", "orange", "blue"), lty=1,bty = "n",cex=0.8 )

dev.off()

latent_infs=c(log10(both_pos$latent_inf),log10(both_neg$latent_inf),log10(hiv_pos$latent_inf),log10(ks_pos$latent_inf))

df<-data.frame(x=latent_infs, y=shedding)

pdf(paste0("arms_",model,"_avgs_latent_inf2.pdf"),height=6,width=6)

ggplot(df, aes(x=x, y=y)) +
  geom_smooth(method=lm,color="black") +
  geom_point(color=arm_colors) +
  xlab("Log reactivation rate") + ylab("Percent Positive") +
  xlim(low_log_latent_inf, high_log_latent_inf) +
  ylim(0, 100) +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12))

df<-data.frame(x=latent_infs, y=log_hhv8, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_latent_inf3.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y)) +
  geom_smooth(method=lm,color="black") +
  geom_point(color=arm_colors) +
  xlab("Log reactivation rate") + ylab("log HHV8 VL") +
  xlim(low_log_latent_inf, high_log_latent_inf) +
  ylim(min(log_hhv8), max(log_hhv8)) +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12))

dev.off()

pdf(paste0("arms_",model,"_avgs_an.pdf"),height=6,width=6)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

# death rate
y_lim=c(low_an,high_an)
ans=list(both_pos$an,both_neg$an,hiv_pos$an,ks_pos$an)

boxplot(ans, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Inf cell death rate (all participants)" ,
	at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
	ylim=y_lim,xlab = "", ylab = "death rate")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

#legend("bottomright", 
#       legend = c("Both Pos", "Both Neg", "HIV Only Pos", "KS Only Pos"), 
#       col = c("red", "green", "orange", "blue"), lty=1,bty = "n" )

ans=c(both_pos$an,both_neg$an,hiv_pos$an,ks_pos$an)

df<-data.frame(x=ans, y=shedding, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_an2.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Percent Positive vs Infected cell death rate") +
  xlab("death rate") + ylab("Percent Positive") +
  xlim(low_an,high_an) +
  ylim(0, 100)

df<-data.frame(x=ans, y=log_hhv8, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_an3.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Median Log HHV8 VL vs Infected cell death rate") +
  xlab("death rate") + ylab("log HHV8 VL") +
  xlim(low_an,high_an) +
  ylim(min(log_hhv8), max(log_hhv8))

dev.off()

pdf(paste0("arms_",model,"_avgs_log_p.pdf"),height=6,width=6)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

#p box plots 
y_lim=c(low_log_p,high_log_p)
p=list((both_pos$log_p),(both_neg$log_p),(hiv_pos$log_p),(ks_pos$log_p))

boxplot(p, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Log Viral Production rate (all participants)" ,
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
        ylim=y_lim,xlab = "", ylab = "log production")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

#legend("bottomright", 
#       legend = c("Both Pos", "Both Neg", "HIV Only Pos", "KS Only Pos"), 
#       col = c("red", "green", "orange", "blue"), lty=1,bty = "n" )

log_p=c((both_pos$log_p),(both_neg$log_p),(hiv_pos$log_p),(ks_pos$log_p))

df<-data.frame(x=log_p, y=shedding, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_log_p2.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Percent Positive vs log Viral Production") +
  xlab("log Viral Production") + ylab("Percent Positive") +
  xlim(low_log_p, high_log_p) +
  ylim(0, 100)

df<-data.frame(x=log_p, y=log_hhv8, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_log_p3.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Median Log HHV8 VL vs log Viral Production") +
  xlab("log Viral Production") + ylab("log HHV8 VL") +
  xlim(low_log_p, high_log_p) +
  ylim(min(log_hhv8), max(log_hhv8))

dev.off()

pdf(paste0("arms_",model,"_avgs_beta.pdf"),height=6,width=6)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

#beta box plots 
y_lim=c(low_beta,high_beta)
beta=list((both_pos$beta),(both_neg$beta),(hiv_pos$beta),(ks_pos$beta))

boxplot(beta, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Cell infectivity (all participants)" ,
	at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
	ylim=y_lim,xlab = "", ylab = "infectivity")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

#legend("bottomright", 
#   legend = c("Both Pos", "Both Neg", "HIV Only Pos", "KS Only Pos"), 
#   col = c("red", "green", "orange", "blue"), lty=1,bty = "n" )

beta=c((both_pos$beta),(both_neg$beta),(hiv_pos$beta),(ks_pos$beta))

df<-data.frame(x=beta, y=shedding, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_beta2.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Percent Positive vs Infectivity") +
  xlab("Infectivity") + ylab("Percent Positive") +
  xlim(low_beta, high_beta) +
  ylim(0, 100)

df<-data.frame(x=beta, y=log_hhv8, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_beta3.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Median Log HHV8 VL vs Infectivity") +
  xlab("Infectivity") + ylab("log HHV8 VL") +
  xlim(low_beta, high_beta) +
  ylim(min(log_hhv8), max(log_hhv8))

dev.off()

pdf(paste0("arms_",model,"_avgs_r.pdf"),height=6,width=6)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

#r box plots 
y_lim=c(low_log_r,high_log_r)
r=list(log10(both_pos$r),log10(both_neg$r),log10(hiv_pos$r),log10(ks_pos$r))

boxplot(r, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Log Tcell Birth IC50 (all participants)" ,
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
        ylim=y_lim,xlab = "", ylab = "log IC50")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

#legend("bottomright", 
#       legend = c("Both Pos", "Both Neg", "HIV Only Pos", "KS Only Pos"), 
#       col = c("red", "green", "orange", "blue"), lty=1,bty = "n" )

log_r=c(log10(both_pos$r),log10(both_neg$r),log10(hiv_pos$r),log10(ks_pos$r))

df<-data.frame(x=log_r, y=shedding, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_r2.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Percent Positive vs Inf Cell IC50") +
  xlab("log IC50") + ylab("Percent Positive") +
  xlim(low_log_r, high_log_r) +
  ylim(0, 100)

df<-data.frame(x=log_r, y=log_hhv8, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_r3.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Median Log HHV8 VL vs Inf Cell IC50 ") +
  xlab("log IC50") + ylab("log HHV8 VL") +
  xlim(low_log_r, high_log_r) +
  ylim(min(log_hhv8), max(log_hhv8))

dev.off()

pdf(paste0("arms_",model,"_avgs_fpos.pdf"),height=6,width=6)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

#fpos box plots 
y_lim=c(low_log_fpos,high_log_fpos)
fpos=list(log10(both_pos$fpos),log10(both_neg$fpos),log10(hiv_pos$fpos),log10(ks_pos$fpos))

boxplot(fpos, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Log Killing rate (all participants)" ,
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
        ylim=y_lim,xlab = "", ylab = "log killing rate")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

#legend("bottomright", 
#       legend = c("Both Pos", "Both Neg", "HIV Only Pos", "KS Only Pos"), 
#       col = c("red", "green", "orange", "blue"), lty=1,bty = "n" )

log_fpos=c(log10(both_pos$fpos),log10(both_neg$fpos),log10(hiv_pos$fpos),log10(ks_pos$fpos))

df<-data.frame(x=log_fpos, y=shedding, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_fpos2.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Percent Positive vs log Killing rate") +
  xlab("log Killing Rate") + ylab("Percent Positive") +
  xlim(low_log_fpos, high_log_fpos) +
  ylim(0, 100)

df<-data.frame(x=log_fpos, y=log_hhv8, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_fpos3.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Median Log HHV8 VL vs log Killing rate") +
  xlab("log killing rate") + ylab("log HHV8 VL") +
  xlim(low_log_fpos, high_log_fpos) +
  ylim(min(log_hhv8), max(log_hhv8))

dev.off()

pdf(paste0("arms_",model,"_avgs_alpha.pdf"),height=6,width=6)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

# birth rate
low_alpha=min(both_pos$alpha,both_neg$alpha,hiv_pos$alpha,ks_pos$alpha)
high_alpha=max(both_pos$alpha,both_neg$alpha,hiv_pos$alpha,ks_pos$alpha)

y_lim=c(low_alpha,high_alpha)
alphas=list(both_pos$alpha,both_neg$alpha,hiv_pos$alpha,ks_pos$alpha)

boxplot(alphas, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Tcell birth rate (all participants)" ,
	at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
	ylim=y_lim,xlab = "", ylab = "birth rate")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

#legend("bottomright", 
#       legend = c("Both Pos", "Both Neg", "HIV Only Pos", "KS Only Pos"), 
#       col = c("red", "green", "orange", "blue"), lty=1,bty = "n" )

alphas=c(both_pos$alpha,both_neg$alpha,hiv_pos$alpha,ks_pos$alpha)

df<-data.frame(x=alphas, y=shedding, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_alpha2.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Percent Positive vs Tcell birth rate") +
  xlab("birth rate") + ylab("Percent Positive") +
  xlim(low_alpha,high_alpha) +
  ylim(0, 100)

df<-data.frame(x=alphas, y=log_hhv8, arm=arms)

dev.off()

pdf(paste0("arms_",model,"_avgs_alpha3.pdf"),height=6,width=6)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  geom_smooth(method=lm) +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  ggtitle("Median Log HHV8 VL vs T Cell birth rate") +
  xlab("birth rate") + ylab("log HHV8 VL") +
  xlim(low_alpha,high_alpha) +
  ylim(min(log_hhv8), max(log_hhv8))

dev.off()
