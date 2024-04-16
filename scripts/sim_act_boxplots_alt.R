library(RColorBrewer)
# Load ggplot2
library(ggplot2)
library(tidyverse)

args<-commandArgs(trailingOnly=T)
ptid_list<-args[1]
desc<-args[2]
model <-args[3]
boxes <-as.numeric(args[4])
leg_pos="none"

height=2.5
width=6.5
if (desc=="both_pos") {
    arm_label="HIV+/KS+"
} else if (desc=="both_neg") {
    arm_label="HIV-/KS-"
} else if (desc=="hiv_pos") {
    arm_label="HIV+/KS-"
} else if (desc=="hiv_pos") {
    arm_label="HIV-/KS+"
} else {
    arm_label="All"
}
ptids  <- file(ptid_list, open = "r")

sim_names=list();
sim_vets=list();
act_names=list();
act_vets=list();

while (length(oneLine <- readLines(ptids, n = 1, warn = FALSE)) > 0) {

    ptid = oneLine;

    sim_data = read.csv(paste0("ptid_",ptid,"_",model,".dat1.csv"))
    act_data = read.csv(paste0("../data/ptid_",ptid,"_act.csv"))
    #time,vet1,vet2,vet3,vet4,vet5,vet6,vet7,vet8,vet9,vet10

    if (sum(sim_data$vet1) > 2) {
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
	print ("No runs with positive swabs (>LOD)! Exiting!")
	quit()
    }
    sim_pos_swabs=list();
    sim_hhv8 = 0
    sim_swabs = 0
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
    sim_hhv8_vl=median(sim_pos_swabs)
    sim_hhv8_shed=100*sim_hhv8 / sim_swabs

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
    act_hhv8_shed=100*act_hhv8 / length(act_data$vet1)

    ptid_num = sub("ptid_", "", ptid)

    if (100*sim_hhv8 / sim_swabs < 20) {
	sim_cat="<20"
    } else if (100*sim_hhv8 / sim_swabs < 40) {
	sim_cat="<40"
    } else if (100*sim_hhv8 / sim_swabs < 60) {
	sim_cat="<60"
    } else if (100*sim_hhv8 / sim_swabs < 80) {
	sim_cat="<80"
    } else {
	sim_cat=">=80"
    }
    if (100*act_hhv8 / length(act_data$vet1) < 20) {
	act_cat="<20"
    } else if (100*act_hhv8 / length(act_data$vet1) < 40) {
	act_cat="<40"
    } else if (100*act_hhv8 / length(act_data$vet1) < 60) {
	act_cat="<60"
    } else if (100*act_hhv8 / length(act_data$vet1) < 80) {
	act_cat="<80"
    } else {
	act_cat=">=80"
    }
    ptid_num = sub("ptid_", "", ptid)
    act_ptid=paste0("A",ptid_num)
    sim_ptid=paste0("S",ptid_num)
    if (length(sim_names)==0) {
	sim_vl_stats = c(rep(act_hhv8_vl,length(sim_pos_swabs)))
	sim_names = c(rep(sim_ptid,length(sim_pos_swabs)))
	sim_groups = c(rep(ptid_num,length(sim_pos_swabs)))
	sim_shed_cats = c(rep(sim_cat,length(sim_pos_swabs)))
	sim_vets = c(sim_pos_swabs)
	act_vl_stats = c(rep(act_hhv8_vl,length(act_pos_swabs)))
	act_names = c(rep(act_ptid,length(act_pos_swabs)))
	act_groups = c(rep(ptid_num,length(act_pos_swabs)))
	act_shed_cats = c(rep(act_cat,length(act_pos_swabs)))
	act_vets = c(act_pos_swabs)
	few_names= c(act_ptid)
	print(paste("first",ptid,"pos swabs",length(sim_pos_swabs)))
    } else {
	sim_vl_stats = c(sim_vl_stats,rep(act_hhv8_vl,length(sim_pos_swabs)))
	sim_names= c(sim_names,rep(sim_ptid,length(sim_pos_swabs)))
	sim_groups = c(sim_groups,rep(ptid_num,length(sim_pos_swabs)))
	sim_shed_cats= c(sim_shed_cats,rep(sim_cat,length(sim_pos_swabs)))
	sim_vets = c(sim_vets,sim_pos_swabs)
	act_vl_stats = c(act_vl_stats,rep(act_hhv8_vl,length(act_pos_swabs)))
	act_names= c(act_names,rep(act_ptid,length(act_pos_swabs)))
	act_groups = c(act_groups,rep(ptid_num,length(act_pos_swabs)))
	act_shed_cats= c(act_shed_cats,rep(act_cat,length(act_pos_swabs)))
	act_vets = c(act_vets,act_pos_swabs)
	few_names= c(few_names,act_ptid)
	print(paste("added",ptid,"pos swabs",length(sim_pos_swabs)))
    }
}
print(paste("Simulated vets:",length(sim_vets),"names:",length(sim_names),"fill:",length(sim_shed_cats),"stats:",length(sim_vl_stats)))
print(paste("checking",length(few_names),"ids vs requested",boxes))
if (length(few_names) < boxes){
    extra_names = seq(1000,1000+(boxes - length(few_names) - 1),1)
    extra_vals = rep(10,(boxes - length(few_names)))
    extra_cats = rep(">80",(boxes - length(few_names)))
    sim_vl_stats = c(sim_vl_stats,extra_vals)
    sim_names= c(sim_names,extra_names)
    sim_groups = c(sim_groups,extra_names)
    sim_shed_cats= c(sim_shed_cats,extra_cats)
    sim_vets = c(sim_vets,extra_vals)
    act_vl_stats = c(act_vl_stats,extra_vals)
    act_names= c(act_names,extra_names)
    act_groups = c(act_groups,extra_names)
    act_shed_cats= c(act_shed_cats,extra_cats)
    act_vets = c(act_vets,extra_vals)
}

y_lim=c(0,8)

print(paste("Simulated vets:",length(sim_vets),"names:",length(sim_names),"fill:",length(sim_shed_cats),"stats:",length(sim_vl_stats)))

pdf(paste0("sim_hhv8_",model,"_",desc,".pdf"),height=height,width=width)

df<-data.frame(vals=sim_vets, names=sim_names, fill=sim_shed_cats, vl_stat=sim_vl_stats, group=sim_groups)

sorted_df <- df[order(df$vl_stat),]
  
colors <- c("<20" = "lightblue", "<40" = "blue", "<60" = "yellow", "<80" = "orange", ">=80" = "red")

ggplot(sorted_df, aes(x=fct_inorder(group),z=fct_inorder(names), y=vals, fill=fill)) +
  #ggtitle("Simulated Shedding by Participant")+
  geom_boxplot() +
  scale_color_manual(values=colors, name="Shedding %",aesthetics = c("fill")) +
  scale_y_continuous(breaks = seq(0, 8, by=1), limits=c(0,8)) +
  xlab("") + ylab("log10 HHV8 DNA") +
  theme(axis.ticks.x = element_blank(),
	   axis.text.x = element_blank(),
	   axis.text.y = element_text(color = "black", size =14),
	   axis.title.y = element_text(color = "black", size =14),
	   legend.position=leg_pos)

dev.off()

y_lim=c(0,8)

print(paste("Actual vets:",length(act_vets),"names:",length(act_names),"fill:",length(act_shed_cats),"stats:",length(act_vl_stats),"group:",length(act_groups)))
pdf(paste0("act_hhv8_",model,"_",desc,".pdf"),height=height,width=width)

df<-data.frame(vals=act_vets, names=act_names, fill=act_shed_cats, vl_stat=act_vl_stats, group=act_groups)

sorted_df <- df[order(df$vl_stat),]
  
colors <- c("<20" = "lightblue", "<40" = "blue", "<60" = "yellow", "<80" = "orange", ">=80" = "red")

ggplot(sorted_df, aes(x=fct_inorder(group),z=fct_inorder(names), y=vals, fill=fill)) +
  #ggtitle("Actual Shedding by Participant")+
  geom_boxplot() +
  scale_color_manual(values=colors, name="Shedding %",aesthetics = c("fill")) +
  scale_y_continuous(breaks = seq(0, 8, by=1), limits=c(0,8)) +
  xlab("") + ylab("log10 HHV8 DNA") +
  theme(axis.ticks.x = element_blank(),
	   axis.text.x = element_blank(),
	   axis.text.y = element_text(color = "black", size =14),
	   axis.title.y = element_text(color = "black", size =14),
	   legend.position="none")

dev.off()
