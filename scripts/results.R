library(RColorBrewer)
args<-commandArgs(trailingOnly=T)
model<-args[1]

both_pos = read.csv(paste0("params_",model,"_both_pos.csv"))
both_neg = read.csv(paste0("params_",model,"_both_neg.csv"))
hiv_pos = read.csv(paste0("params_",model,"_hiv_pos.csv"))
ks_pos = read.csv(paste0("params_",model,"_ks_pos.csv"))
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
high_log_latent_inf=2

low_exp_days=7;
high_exp_days=20;

low_log_r=1;
high_log_r=3;

if (model==5) {
    pdf(paste0("betaI_model_",model,"_misc_params.pdf"),height=6,width=12)
    par(mfrow = c(1,3), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))
} else {
    pdf(paste0("betaI_model_",model,"_misc_params.pdf"),height=6,width=6)
    par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))
}

# latent_reacts
y_lim=c(low_log_latent_inf,high_log_latent_inf)
latent_infs=list(log10(both_pos$latent_inf),log10(both_neg$latent_inf),log10(hiv_pos$latent_inf),log10(ks_pos$latent_inf))

boxplot(latent_infs, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Log Activation rate",
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
        ylim=y_lim,xlab = "", ylab = "log latent inf")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

legend("bottomright", 
       legend = c("Both Pos", "Both Neg", "HIV Only Pos", "KS Only Pos"), 
       col = c("red", "green", "orange", "blue"), lty=1,bty = "n" )

if (model==5) {
# death rate
y_lim=c(low_an,high_an)
ans=list(both_pos$an,both_neg$an,hiv_pos$an,ks_pos$an)

boxplot(ans, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Inf cell death rate",
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
        ylim=y_lim,xlab = "", ylab = "an")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

legend("bottomright", 
       legend = c("Both Pos", "Both Neg", "HIV Only Pos", "KS Only Pos"), 
       col = c("red", "green", "orange", "blue"), lty=1,bty = "n" )
}

dev.off()

if (model==5) {
    pdf(paste0("betaI_model_",model,"_viral_params.pdf"),height=6,width=8)
    par(mfrow = c(1,2), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))
} else {
    pdf(paste0("betaI_model_",model,"_viral_params.pdf"),height=6,width=6)
    par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))
}

if (model==5) {
#beta box plots 
y_lim=c(low_beta,high_beta)
beta=list((both_pos$beta),(both_neg$beta),(hiv_pos$beta),(ks_pos$beta))

boxplot(beta, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="cell infectivity",
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
        ylim=y_lim,xlab = "", ylab = "log beta")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)
}

#p box plots 
y_lim=c(low_log_p,high_log_p)
p=list((both_pos$log_p),(both_neg$log_p),(hiv_pos$log_p),(ks_pos$log_p))

boxplot(p, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Log Production rate",
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
        ylim=y_lim,xlab = "", ylab = "log p")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

legend("bottomright", 
       legend = c("Both Pos", "Both Neg", "HIV Only Pos", "KS Only Pos"), 
       col = c("red", "green", "orange", "blue"), lty=1,bty = "n" )

dev.off()

pdf(paste0("betaI_model_",model,"_tcell_params.pdf"),height=6,width=8)
par(mfrow = c(1,2), mar = 0.1 + c(3, 3, 1, 1), mgp = c(2, 0.5, 0))

#r box plots 
y_lim=c(low_log_r,high_log_r)
r=list(log10(both_pos$r),log10(both_neg$r),log10(hiv_pos$r),log10(ks_pos$r))

boxplot(r, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Log Tcell Birth IC50",
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
        ylim=y_lim,xlab = "", ylab = "log r")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)


#fpos box plots 
y_lim=c(low_log_fpos,high_log_fpos)
fpos=list(log10(both_pos$fpos),log10(both_neg$fpos),log10(hiv_pos$fpos),log10(ks_pos$fpos))

boxplot(fpos, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Log Killing rate",
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim = c(0, 7.0),
        ylim=y_lim,xlab = "", ylab = "log fpos")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

legend("bottomright", 
       legend = c("Both Pos", "Both Neg", "HIV Only Pos", "KS Only Pos"), 
       col = c("red", "green", "orange", "blue"), lty=1,bty = "n" )

dev.off()
