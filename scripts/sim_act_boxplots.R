library(RColorBrewer)
# Load ggplot2
library(ggplot2)
library(tidyverse)

args<-commandArgs(trailingOnly=T)
ptid_list<-args[1]
desc<-args[2]
model <-args[3]

ptids  <- file(ptid_list, open = "r")

names=list();
vets=list();
alphas=list();

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
	sim_alpha="Simulated"
	sim_cat="<20"
    } else if (100*sim_hhv8 / sim_swabs < 40) {
	sim_alpha="Simulated"
	sim_cat="<40"
    } else if (100*sim_hhv8 / sim_swabs < 60) {
	sim_alpha="Simulated"
	sim_cat="<60"
    } else if (100*sim_hhv8 / sim_swabs < 80) {
	sim_alpha="Simulated"
	sim_cat="<80"
    } else {
	sim_alpha="Simulated"
	sim_cat=">=80"
    }
    if (100*act_hhv8 / length(act_data$vet1) < 20) {
	act_alpha="Actual"
	act_cat="<20"
    } else if (100*act_hhv8 / length(act_data$vet1) < 40) {
	act_alpha="Actual"
	act_cat="<40"
    } else if (100*act_hhv8 / length(act_data$vet1) < 60) {
	act_alpha="Actual"
	act_cat="<60"
    } else if (100*act_hhv8 / length(act_data$vet1) < 80) {
	act_alpha="Actual"
	act_cat="<80"
    } else {
	act_alpha="Actual"
	act_cat=">=80"
    }
    ptid_num = sub("ptid_", "", ptid)
    act_ptid=paste0("A",ptid_num)
    sim_ptid=paste0("S",ptid_num)
    if (length(names)==0) {
	vl_stats = c(rep(act_hhv8_shed,length(sim_pos_swabs)),rep(act_hhv8_shed,length(act_pos_swabs)))
	few_vl_stats = c(act_hhv8_shed,act_hhv8_shed)
	shed_stats = c(sim_hhv8_shed,act_hhv8_shed)
	few_names= c(sim_ptid,act_ptid)
	few_groups= c(ptid_num,ptid_num)
	names = c(rep(sim_ptid,length(sim_pos_swabs)),rep(act_ptid,length(act_pos_swabs)))
	groups = c(rep(ptid_num,length(sim_pos_swabs)),rep(ptid_num,length(act_pos_swabs)))
	alphas = c(rep(sim_alpha,length(sim_pos_swabs)),rep(act_alpha,length(act_pos_swabs)))
	few_alphas = c(sim_alpha,act_alpha)
	shed_cats = c(rep(sim_cat,length(sim_pos_swabs)),rep(act_cat,length(act_pos_swabs)))
	few_shed_cats = c(sim_cat,act_cat)
	vets = c(sim_pos_swabs,act_pos_swabs)
	print(paste("first",ptid))
    } else {
	vl_stats = c(vl_stats,rep(act_hhv8_shed,length(sim_pos_swabs)),rep(act_hhv8_shed,length(act_pos_swabs)))
	few_vl_stats = c(few_vl_stats,act_hhv8_shed,act_hhv8_shed)
	shed_stats = c(shed_stats,sim_hhv8_shed,act_hhv8_shed)
	few_names= c(few_names,sim_ptid,act_ptid)
	names= c(names,rep(sim_ptid,length(sim_pos_swabs)),rep(act_ptid,length(act_pos_swabs)))
	groups = c(groups,rep(ptid_num,length(sim_pos_swabs)),rep(ptid_num,length(act_pos_swabs)))
	few_groups= c(few_groups,ptid_num,ptid_num)
	alphas = c(alphas,rep(sim_alpha,length(sim_pos_swabs)),rep(act_alpha,length(act_pos_swabs)))
	few_alphas = c(few_alphas,sim_alpha,act_alpha)
	shed_cats= c(shed_cats,rep(sim_cat,length(sim_pos_swabs)),rep(act_cat,length(act_pos_swabs)))
	few_shed_cats = c(few_shed_cats,sim_cat,act_cat)
	vets = c(vets,sim_pos_swabs,act_pos_swabs)
	print(paste("added",ptid))
    }
}

y_lim=c(0,8)

pdf(paste0("sim_vs_act_hhv8_",model,"_",desc,".pdf"),height=8,width=16)

df<-data.frame(vals=vets, names=names, fill=shed_cats, vl_stat=vl_stats, Data=alphas,group=groups)

sorted_df <- df[order(df$vl_stat),]
options(max.print=999999)
sorted_df
  
colors <- c("<20" = "lightblue", "<40" = "blue", "<60" = "yellow", "<80" = "orange", ">=80" = "red")

ggplot(sorted_df, aes(x=fct_inorder(group),z=fct_inorder(names), y=vals, fill=fill, alpha=Data)) +
  geom_boxplot() +
  #ggtitle(paste0("Log HHV8 VL Sim & Act (",desc," model ",model,")")) +
  scale_color_manual(values=colors, name="Shedding %",aesthetics = c("fill")) +
  scale_alpha_discrete(range = c(0.5, 1.0)) +
  xlab("PTID") + ylab("log HHV8")+
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 14))

dev.off()

pdf(paste0("sim_vs_act_shedding_",model,"_",desc,".pdf"),height=8,width=16)

df<-data.frame(vals=shed_stats, names=few_names, fill=few_shed_cats, vl_stat=few_vl_stats, Data=few_alphas,group=few_groups)

sorted_df <- df[order(df$vl_stat),]

ggplot(sorted_df, aes(x=fct_inorder(group), z=fct_inorder(names), y=vals, fill=fill, alpha=Data)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  #ggtitle(paste0("Shedding Percent Sim & Act (",desc," model ",model,")")) +
  scale_color_manual(values=colors, name="Shedding %",aesthetics = c("fill")) +
  scale_alpha_discrete(range = c(0.5, 1.0)) +
  xlab("PTID") + ylab("Percent Shedding")+
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 14))

dev.off()
