library(RColorBrewer)
# Load ggplot2
library(ggplot2)
library(tidyverse)
library(viridis)


args<-commandArgs(trailingOnly=T)
model<-args[1]
desc<-args[2]

stats_hiv = read.csv(paste0("act_",model,"_all_hiv.csv"))
stats_no_hiv = read.csv(paste0("act_",model,"_no_hiv.csv"))
stats_ks = read.csv(paste0("act_",model,"_all_ks.csv"))
stats_no_ks = read.csv(paste0("act_",model,"_no_ks.csv"))

stats_both_pos = read.csv(paste0("act_",model,"_both_pos.csv"))
stats_both_neg = read.csv(paste0("act_",model,"_both_neg.csv"))
stats_hiv_pos = read.csv(paste0("act_",model,"_hiv_pos.csv"))
stats_ks_pos = read.csv(paste0("act_",model,"_ks_pos.csv"))

hiv_shedding=c(stats_hiv$shedding,stats_no_hiv$shedding)
ks_shedding=c(stats_ks$shedding,stats_no_ks$shedding)
hiv_med_hhv8=c(stats_hiv$med_hhv8,stats_no_hiv$med_hhv8)
ks_med_hhv8=c(stats_ks$med_hhv8,stats_no_ks$med_hhv8)
hiv_peak_hhv8=c(stats_hiv$peak_hhv8,stats_no_hiv$peak_hhv8)
ks_peak_hhv8=c(stats_ks$peak_hhv8,stats_no_ks$peak_hhv8)

hiv_names=c(rep("HIV Pos",length(stats_hiv$shedding)),rep("HIV Neg",length(stats_no_hiv$shedding)))
ks_names=c(rep("KS Pos",length(stats_ks$shedding)),rep("KS Neg",length(stats_no_ks$shedding)))
arm_names=c(rep("Both Pos",length(stats_both_pos$shedding)),rep("Both Neg",length(stats_both_neg$shedding)),rep("HIV Pos",length(stats_hiv_pos$shedding)),rep("KS Pos",length(stats_ks_pos$shedding)))

pdf(paste0(desc,"_shed_box_HIV.pdf"),height=3,width=3)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 3, 1), mgp = c(2, 0.5, 0))

df<-data.frame(x=hiv_names, y=hiv_shedding,Status=hiv_names)
ggplot(df, aes(x=x, y=y, color=Status)) +
  geom_boxplot() +
  scale_color_manual(values = c("HIV Pos" = "red", "HIV Neg" = "green")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("Shedding by HIV Status") +
  xlab("") + ylab("Percent Positive") +
  ylim(0, 100) + theme(legend.position = "none")
dev.off()

pdf(paste0(desc,"_shed_box_KS.pdf"),height=3,width=3)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 3, 1), mgp = c(2, 0.5, 0))

df<-data.frame(x=ks_names, y=ks_shedding,Status=ks_names)
ggplot(df, aes(x=x, y=y, color=Status)) +
  geom_boxplot() +
  scale_color_manual(values = c("KS Pos" = "red", "KS Neg" = "green")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("Shedding by KS Status") +
  xlab("") + ylab("Percent Positive") +
  ylim(0, 100) + theme(legend.position = "none")
dev.off()

#boxplot(hiv_shedding, range=0.0, horizontal=FALSE, col=c("red", "green"),
#	main="Percent Positive Shedding by HIV Status",
#        at=c(2,4),  beside=TRUE, xaxt = "n",xlim=c(0,6),
#        ylim=c(0, 100.0),xlab = "", ylab = "percent positive")
#axis(side = 1, line = 0, at = c(2, 4), labels = c("HIV Pos", "HIV Neg"), tick = FALSE, cex.axis = 0.8)
#
#boxplot(ks_shedding, range=0.0, horizontal=FALSE, col=c("red", "green"),
#	main="Percent Positive Shedding by KS Status",
#        at=c(2,4),  beside=TRUE, xaxt = "n",xlim=c(0,6),
#        ylim=c(0, 100.0),xlab = "", ylab = "percent positive")
#axis(side = 1, line = 0, at = c(2, 4), labels = c("KS Pos", "KS Neg"), tick = FALSE, cex.axis = 0.8)

pdf(paste0(desc,"_shed_box_ARMs.pdf"),height=3,width=3)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 3, 1), mgp = c(2, 0.5, 0))

shedding=c(stats_both_pos$shedding,stats_both_neg$shedding,stats_hiv_pos$shedding,stats_ks_pos$shedding)
df<-data.frame(x=arm_names, y=shedding,Status=arm_names)
ggplot(df, aes(x=x, y=y, color=Status)) +
  geom_boxplot() +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange","KS Pos"="blue")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("Shedding by Study ARM") +
  xlab("") + ylab("Percent Positive") +
  ylim(0, 100) + theme(legend.position = "none")
#shedding=list(stats_both_pos$shedding,stats_both_neg$shedding,stats_hiv_pos$shedding,stats_ks_pos$shedding)
#boxplot(shedding, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
#	main="Percent Positive Shedding by ARM\n(Participants w/ Pos Swab)",
#        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim=c(0,7),
#        ylim=c(0, 100.0),xlab = "", ylab = "percent positive")
#axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

dev.off()

pdf(paste0(desc,"_shed_by_HIV.pdf"),height=3,width=3)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 3, 1), mgp = c(2, 0.5, 0))
hiv_shedding=c(sum(stats_hiv$shedding)/length(stats_hiv$shedding),sum(stats_no_hiv$shedding)/length(stats_no_hiv$shedding))
hiv_shedding

barplot(hiv_shedding, beside=TRUE, col=c("red", "green"),
	main="Percent Positive Shedding by HIV Status\n(Participants w/ Pos Swab)",
	names.arg=c("HIV Pos", "HIV Neg"),
        #at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim=c(0,7),
        ylim=c(0, 100.0),xlab = "", ylab = "percent positive")
#axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)
dev.off()

pdf(paste0(desc,"_shed_by_KS.pdf"),height=3,width=3)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 3, 1), mgp = c(2, 0.5, 0))
ks_shedding=c(sum(stats_ks$shedding)/length(stats_ks$shedding),sum(stats_no_ks$shedding)/length(stats_no_ks$shedding))
ks_shedding

barplot(ks_shedding, beside=TRUE, col=c("red", "green"),
	main="Percent Positive Shedding by KS Status\n(Participants w/ Pos Swab)",
	names.arg=c("KS Pos", "KS Neg"),
        #at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim=c(0,7),
        ylim=c(0, 100.0),xlab = "", ylab = "percent positive")
#axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)
dev.off()

pdf(paste0(desc,"_shed_by_ARM.pdf"),height=3,width=3)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 3, 1), mgp = c(2, 0.5, 0))
shedding=c(sum(stats_both_pos$shedding)/length(stats_both_pos$shedding),sum(stats_both_neg$shedding)/length(stats_both_neg$shedding),sum(stats_hiv_pos$shedding)/length(stats_hiv_pos$shedding),sum(stats_ks_pos$shedding)/length(stats_ks_pos$shedding))
shedding

barplot(shedding, beside=TRUE, col=c("red", "green", "orange", "blue"),
	main="Percent Positive Shedding by ARM\n(Participants w/ Pos Swab)",
	names.arg=c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"),
        #at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim=c(0,7),
        ylim=c(0, 100.0),xlab = "", ylab = "percent positive")
#axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

dev.off()

hiv_names=c(rep("HIV Pos",length(stats_hiv$shedding)),rep("HIV Neg",length(stats_no_hiv$shedding)))
ks_names=c(rep("KS Pos",length(stats_ks$shedding)),rep("KS Neg",length(stats_no_ks$shedding)))


pdf(paste0(desc,"_med_hhv8_box_HIV.pdf"),height=3,width=3)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 3, 1), mgp = c(2, 0.5, 0))

df<-data.frame(x=hiv_names, y=hiv_med_hhv8,Status=hiv_names)
ggplot(df, aes(x=x, y=y, color=Status)) +
  geom_boxplot() +
  scale_color_manual(values = c("HIV Pos" = "red", "HIV Neg" = "green")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("Median log HHV8 VL vs HIV") +
  xlab("") + ylab("log HHV8") +
  ylim(0, 8) + theme(legend.position = "none")

dev.off()

pdf(paste0(desc,"_med_hhv8_box_KS.pdf"),height=3,width=3)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 3, 1), mgp = c(2, 0.5, 0))

df<-data.frame(x=ks_names, y=ks_med_hhv8,Status=ks_names)
ggplot(df, aes(x=x, y=y, color=Status)) +
  geom_boxplot() +
  scale_color_manual(values = c("KS Pos" = "red", "KS Neg" = "green")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("Median log HHV8 VL vs KS") +
  xlab("") + ylab("log HHV8") +
  ylim(0, 8) + theme(legend.position = "none")

med_hhv8=c(stats_both_pos$med_hhv8,stats_both_neg$med_hhv8,stats_hiv_pos$med_hhv8,stats_ks_pos$med_hhv8)

dev.off()

pdf(paste0(desc,"_med_hhv8_box_ARMs.pdf"),height=3,width=3)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 3, 1), mgp = c(2, 0.5, 0))

df<-data.frame(x=arm_names, y=med_hhv8,Status=arm_names)
ggplot(df, aes(x=x, y=y, color=Status)) +
  geom_boxplot() +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange","KS Pos"="blue")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("Median log HHV8 VL vs ARM") +
  xlab("") + ylab("log HHV8") +
  ylim(0, 8) + theme(legend.position = "none")

#boxplot(med_hhv8, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
#	main="Median log HHV8 VL by ARM\n(Participants w/ Pos Swab)",
#        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim=c(0,7),
#        ylim=c(0, 8.0),xlab = "", ylab = "log VL")
#axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

dev.off()

pdf(paste0(desc,"_peak_hhv8_vl.pdf"),height=3,width=3)
par(mfrow = c(1,1), mar = 0.1 + c(3, 3, 3, 1), mgp = c(2, 0.5, 0))

df<-data.frame(x=hiv_names, y=hiv_peak_hhv8,Status=hiv_names)
ggplot(df, aes(x=x, y=y, color=Status)) +
  geom_boxplot() +
  scale_color_manual(values = c("HIV Pos" = "red", "HIV Neg" = "green")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("Peak log HHV8 VL by HIV Status") +
  xlab("") + ylab("log HHV8") +
  ylim(0, 8) + theme(legend.position = "none")


df<-data.frame(x=ks_names, y=ks_peak_hhv8,Status=ks_names)
ggplot(df, aes(x=x, y=y, color=Status)) +
  geom_boxplot() +
  scale_color_manual(values = c("KS Pos" = "red", "KS Neg" = "green")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle("Peak log HHV8 VL by KS Status") +
  xlab("") + ylab("log HHV8") +
  ylim(0, 8) + theme(legend.position = "none")

peak_hhv8=list(stats_both_pos$peak_hhv8,stats_both_neg$peak_hhv8,stats_hiv_pos$peak_hhv8,stats_ks_pos$peak_hhv8)

boxplot(peak_hhv8, range=0.0, horizontal=FALSE, col=c("red", "green", "orange", "blue"),
	main="Peak log HHV8 VL by ARM\n(Participants w/ Pos Swab)",
        at=c(1,2.5,4,5.5),  beside=TRUE, xaxt = "n",xlim=c(0,7),
        ylim=c(0, 8.0),xlab = "", ylab = "log VL")
axis(side = 1, line = 0, at = c(1, 2.5, 4, 5.5), labels = c("Both Pos", "Both Neg", "HIV Pos", "KS Pos"), tick = FALSE, cex.axis = 0.8)

dev.off()
