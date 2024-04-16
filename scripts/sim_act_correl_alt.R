library(RColorBrewer)
# Load ggplot2
library(ggplot2)
library(tidyverse)

args<-commandArgs(trailingOnly=T)
desc <-args[1]
model <-args[2]
both_pos_file = file("../data/both_pos_ptids.txt",open="r")
both_neg_file = file("../data/both_neg_ptids.txt",open="r")
hiv_pos_file = file("../data/hiv_pos_ptids.txt",open="r")
ks_pos_file = file("../data/ks_pos_ptids.txt",open="r")

both_pos=list()
i=1
while (length(oneLine <- readLines(both_pos_file, n = 1, warn = FALSE)) > 0) {
    both_pos[i] = oneLine;
    i=i+1
}

both_neg=list()
i=1
while (length(oneLine <- readLines(both_neg_file, n = 1, warn = FALSE)) > 0) {
    both_neg[i] = oneLine;
    i=i+1
}

hiv_pos=list()
i=1
while (length(oneLine <- readLines(hiv_pos_file, n = 1, warn = FALSE)) > 0) {
    hiv_pos[i] = oneLine;
    i=i+1
}

ks_pos=list()
i=1
while (length(oneLine <- readLines(ks_pos_file, n = 1, warn = FALSE)) > 0) {
    ks_pos[i] = oneLine;
    i=i+1
}

ptids  <- c(both_pos,both_neg,hiv_pos,ks_pos)

arms=c(rep("Both Pos",length(both_pos)),rep("Both Neg",length(both_neg)),rep("HIV Pos",length(hiv_pos)),rep("KS Pos",length(ks_pos)))

sim_shed_perc=list();
act_shed_perc=list();
sim_med_hhv8=list();
act_med_hhv8=list();
sim_max_hhv8=list();
act_max_hhv8=list();

first_line=1;

write("ptid,shedding,med_hhv8,peak_hhv8",file=paste0("sim_",model,"_",desc,".csv"),append=F)
write("ptid,shedding,med_hhv8,peak_hhv8",file=paste0("act_",model,"_",desc,".csv"),append=F)
for (i in 1:length(ptids)) {

    ptid = ptids[i];

    sim_data = read.csv(paste0("ptid_",ptid,"_",model,".dat1.csv"))
    #time,vet1,vet2,vet3,vet4,vet5,vet6,vet7,vet8,vet9,vet10

    if (max(sim_data$vet1) > 2) {
	sim_vet = sim_data$vet1
    } else if (max(sim_data$vet2) > 2) {
	sim_vet = sim_data$vet2
    } else if (max(sim_data$vet3) > 2) {
	sim_vet = sim_data$vet3
    } else if (max(sim_data$vet4) > 2) {
	sim_vet = sim_data$vet4
    } else if (max(sim_data$vet5) > 2) {
	sim_vet = sim_data$vet5
    } else if (max(sim_data$vet6) > 2) {
	sim_vet = sim_data$vet6
    } else if (max(sim_data$vet7) > 2) {
	sim_vet = sim_data$vet7
    } else if (max(sim_data$vet8) > 2) {
	sim_vet = sim_data$vet8
    } else if (max(sim_data$vet9) > 2) {
	sim_vet = sim_data$vet9
    } else if (max(sim_data$vet10) > 2) {
	sim_vet = sim_data$vet10
    } else {
	print (paste0("No runs with positive swabs for ptid_",ptid,"_",model,"! Exiting!"))
	quit()
    }
    first_qual=0
    ptid_shed_perc=list();
    ptid_med_hhv8s=list();
    ptid_max_hhv8s=list();

    for (j in 1:10) {

	if (j==1 && max(sim_data$vet1) > 2) {
	    sim_vet = sim_data$vet1
	} else if (j==2 && max(sim_data$vet2) > 2) {
	    sim_vet = sim_data$vet2
	} else if (j==3 && max(sim_data$vet3) > 2) {
	    sim_vet = sim_data$vet3
	} else if (j==4 && max(sim_data$vet4) > 2) {
	    sim_vet = sim_data$vet4
	} else if (j==5 && max(sim_data$vet5) > 2) {
	    sim_vet = sim_data$vet5
	} else if (j==6 && max(sim_data$vet6) > 2) {
	    sim_vet = sim_data$vet6
	} else if (j==7 && max(sim_data$vet7) > 2) {
	    sim_vet = sim_data$vet7
	} else if (j==8 && max(sim_data$vet8) > 2) {
	    sim_vet = sim_data$vet8
	} else if (j==9 && max(sim_data$vet9) > 2) {
	    sim_vet = sim_data$vet9
	} else if (j==10 && max(sim_data$vet10) > 2) {
	    sim_vet = sim_data$vet10
	} else {
	    next;
	} 
	sim_hhv8 = 0
	sim_swabs = 0
	sim_pos_swabs=list();
	offset_first_pos=0.5
	# downsample!
	# find first detectable swab, then swab one day later
	for (i in 1:length(sim_vet)) {
	    if (sim_vet[i] > 2) {
		print (paste("1st positive swabs (>LOD) at",sim_data$time[i]))
		offset_first_pos = sim_data$time[i]-as.integer(sim_data$time[i])
		break;
	    }
	}
	for (i in 1:length(sim_vet)) {
	    if ((sim_data$time[i]-as.integer(sim_data$time[i]) > (offset_first_pos-0.01)) && 
		    (sim_data$time[i]-as.integer(sim_data$time[i]) < (offset_first_pos+0.01))) {
		if (sim_vet[i] > 2) {
		    sim_hhv8 = sim_hhv8+1;
		    sim_pos_swabs= c(sim_pos_swabs,sim_vet[i])
		}
		sim_swabs=sim_swabs+1;
	    }
	}
	sim_pos_swabs <- unlist(sim_pos_swabs)
	sim_pos_swabs <- as.vector(sim_pos_swabs,'numeric')

	this_hhv8_vl=median(sim_pos_swabs)
	this_hhv8_max=max(sim_pos_swabs)

	this_hhv8_shed=100*sim_hhv8 / sim_swabs
	print (paste(this_hhv8_shed,"perc positive swabs in run",j))

	if (first_qual==0) {
	    ptid_shed_percs = this_hhv8_shed
	    ptid_med_hhv8s = this_hhv8_vl
	    ptid_max_hhv8s = this_hhv8_max
	    first_qual=1
	} else {
	    ptid_shed_percs= c(ptid_shed_percs,this_hhv8_shed)
	    ptid_med_hhv8s= c(ptid_med_hhv8s,this_hhv8_vl)
	    ptid_max_hhv8s= c(ptid_max_hhv8s,this_hhv8_max)
	}
    }
    sim_hhv8_shed = mean(ptid_shed_percs)
    sim_hhv8_vl = mean(ptid_med_hhv8s)
    sim_hhv8_max = mean(ptid_max_hhv8s)
    print (paste(sim_hhv8_shed,"perc positive swabs, med log HHV8 =",sim_hhv8_vl,"max =",sim_hhv8_max))

    write(paste(ptid,sim_hhv8_shed,sim_hhv8_vl,sep=","),file=paste0("sim_",model,"_",desc,".csv"),append=T)
    act_data = read.csv(paste0("./ptid_",ptid,"_act.csv"))
    act_hhv8 = 0
    act_pos_swabs=list();
    for (i in 1:length(act_data$vet1)) {
	if (act_data$vet1[i] > 2) {
	    act_hhv8 = act_hhv8+1;
	    act_pos_swabs= c(act_pos_swabs,act_data$vet1[i])
	}
    }
    act_pos_swabs <- unlist(act_pos_swabs)
    act_pos_swabs <- as.vector(act_pos_swabs,'numeric')
    act_hhv8_vl=median(act_pos_swabs)
    act_hhv8_max=max(act_pos_swabs)
    act_hhv8_shed=100*act_hhv8 / length(act_data$vet1)
    print (paste(act_hhv8_shed,"perc positive swabs, med log HHV8 =",act_hhv8_vl,"max =",act_hhv8_max))
    write(paste(ptid,act_hhv8_shed,act_hhv8_vl,sep=","),file=paste0("act_",model,"_",desc,".csv"),append=T)

    if (first_line==1) {
	sim_shed_perc = sim_hhv8_shed
	act_shed_perc = act_hhv8_shed
	sim_med_hhv8 = sim_hhv8_vl
	sim_max_hhv8 = sim_hhv8_max
	act_med_hhv8 = act_hhv8_vl
	act_max_hhv8 = act_hhv8_max
    } else {
	sim_shed_perc= c(sim_shed_perc,sim_hhv8_shed)
	act_shed_perc= c(act_shed_perc,act_hhv8_shed)
	sim_med_hhv8= c(sim_med_hhv8,sim_hhv8_vl)
	act_med_hhv8= c(act_med_hhv8,act_hhv8_vl)
	sim_max_hhv8= c(sim_max_hhv8,sim_hhv8_max)
	act_max_hhv8= c(act_max_hhv8,act_hhv8_max)
    }
    print(paste("added",ptid))
    first_line=0
}

pdf(paste0("model_",model,"_shed_perc_",desc,".pdf"),height=8,width=8)

print("Calculating concordance correlation coefficient for perc positive swabs")
rho=cor(act_shed_perc, sim_shed_perc)
mean1=mean(act_shed_perc)
std1=sd(act_shed_perc)
mean2=mean(sim_shed_perc)
std2=sd(sim_shed_perc)
concord = 2*rho*std1*std2/(std1*std1 + std2*std2 + (mean1-mean2)**2)
concord

df<-data.frame(x=act_shed_perc, y=sim_shed_perc,arm=arms)

lower=min(act_shed_perc,sim_shed_perc)
upper=max(act_shed_perc,sim_shed_perc)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  geom_smooth(method=lm) +
  geom_abline(intercept = 0, slope = 1,color="red", 
                 linetype="dashed")+
  ggtitle(paste0("Shedding Percent Sim vs Actual (model ",model," ",desc,")")) +
  xlab("Actual") + ylab("Simulated") +
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 14)) +
  xlim(lower, upper) +
  ylim(lower, upper)

dev.off()

pdf(paste0("model_",model,"_hhv8_vl_",desc,".pdf"),height=8,width=8)

print("Calculating concordance correlation coefficient for median HHV8 swab")
rho=cor(act_med_hhv8, sim_med_hhv8)
mean1=mean(act_med_hhv8)
std1=sd(act_med_hhv8)
mean2=mean(sim_med_hhv8)
std2=sd(sim_med_hhv8)
concord = 2*rho*std1*std2/(std1*std1 + std2*std2 + (mean1-mean2)**2)
concord

df<-data.frame(x=act_med_hhv8, y=sim_med_hhv8, arm=arms)

lower=min(act_med_hhv8,sim_med_hhv8)
upper=max(act_med_hhv8,sim_med_hhv8)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  geom_smooth(method=lm) +
  geom_abline(intercept = 0, slope = 1,color="red", 
                 linetype="dashed")+
  ggtitle(paste0("Median Log HHV8 VL Sim & Actual (model ",model," ",desc,")")) +
  xlab("Actual") + ylab("Simulated") +
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 14)) +
  xlim(lower, upper) +
  ylim(lower, upper)

dev.off()

pdf(paste0("model_",model,"_hhv8_peak_",desc,".pdf"),height=8,width=8)

print("Calculating concordance correlation coefficient for peak HHV8 swab")
rho=cor(act_max_hhv8, sim_max_hhv8)
mean1=mean(act_max_hhv8)
std1=sd(act_max_hhv8)
mean2=mean(sim_max_hhv8)
std2=sd(sim_max_hhv8)
concord = 2*rho*std1*std2/(std1*std1 + std2*std2 + (mean1-mean2)**2)
concord

df<-data.frame(x=act_max_hhv8, y=sim_max_hhv8, arm=arms)

lower=min(act_max_hhv8,sim_max_hhv8)
upper=max(act_max_hhv8,sim_max_hhv8)
ggplot(df, aes(x=x, y=y, color=arm)) +
  geom_point() +
  scale_color_manual(values = c("Both Pos" = "red", "Both Neg" = "green", "HIV Pos" = "orange", "KS Pos" = "blue")) +
  geom_smooth(method=lm) +
  geom_abline(intercept = 0, slope = 1,color="red", 
                 linetype="dashed")+
  ggtitle(paste0("Peak Log HHV8 VL Sim & Actual (model ",model," ",desc,")")) +
  xlab("Actual") + ylab("Simulated") +
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 14)) +
  xlim(lower, upper) +
  ylim(lower, upper)

dev.off()
