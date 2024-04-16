library(RColorBrewer)
# Load ggplot2
library(ggplot2)
library(tidyverse)

args<-commandArgs(trailingOnly=T)
ptid_list<-args[1]
desc<-args[2]
model <-args[3]

ptids  <- file(ptid_list, open = "r")
list_root<- substr(ptid_list,1,nchar(ptid_list)-10)

names=list();
vets=list();

while (length(oneLine <- readLines(ptids, n = 1, warn = FALSE)) > 0) {

    ptid = oneLine;

    act_data = read.csv(paste0("./ptid_",ptid,"_act.csv"))
    #time,vet1,vet2,vet3,vet4,vet5,vet6,vet7,vet8,vet9,vet10

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
    if (length(names)==0) {
	vl_stats = c(rep(act_hhv8_vl,length(act_pos_swabs)))
	few_vl_stats = c(act_hhv8_vl)
	shed_stats = c(act_hhv8_shed)
	few_names= c(act_ptid)
	few_groups= c(ptid_num)
	names = c(rep(act_ptid,length(act_pos_swabs)))
	groups = c(rep(ptid_num,length(act_pos_swabs)))
	shed_cats = c(rep(act_cat,length(act_pos_swabs)))
	few_shed_cats = c(act_cat)
	vets = c(act_pos_swabs)
	print(paste("first",ptid))
    } else {
	vl_stats = c(vl_stats,rep(act_hhv8_vl,length(act_pos_swabs)))
	few_vl_stats = c(few_vl_stats,act_hhv8_vl)
	shed_stats = c(shed_stats,act_hhv8_shed)
	few_names= c(few_names,act_ptid)
	names= c(names,rep(act_ptid,length(act_pos_swabs)))
	groups = c(groups,rep(ptid_num,length(act_pos_swabs)))
	few_groups= c(few_groups,ptid_num)
	shed_cats= c(shed_cats,rep(act_cat,length(act_pos_swabs)))
	few_shed_cats = c(few_shed_cats,act_cat)
	vets = c(vets,act_pos_swabs)
	print(paste("added",ptid))
    }
}

y_lim=c(0,8)

pdf(paste0("act_hhv8_",model,"_",list_root,".pdf"),height=2.5,width=11)

df<-data.frame(vals=vets, names=names, fill=shed_cats, vl_stat=vl_stats, group=groups)

sorted_df <- df[order(df$vl_stat),]
options(max.print=999999)
sorted_df
  
colors <- c("<20" = "lightblue", "<40" = "blue", "<60" = "yellow", "<80" = "orange", ">=80" = "red")

ggplot(sorted_df, aes(x=fct_inorder(group),z=fct_inorder(names), y=vals, fill=fill)) +
  geom_boxplot() +
  ggtitle(paste0("Log HHV8 VL Swabs for ",desc," Study ARM")) +
  scale_color_manual(values=colors, name="Shedding %",aesthetics = c("fill")) +
  ylim(0,8) +
  xlab("PTIDs") + ylab("log HHV8 VL")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

dev.off()

pdf(paste0("act_shedding_",model,"_",list_root,".pdf"),height=2.5,width=8)

df<-data.frame(vals=shed_stats, names=few_names, fill=few_shed_cats, vl_stat=few_vl_stats, group=few_groups)

sorted_df <- df[order(df$vl_stat),]

ggplot(sorted_df, aes(x=fct_inorder(group), z=fct_inorder(names), y=vals, fill=fill)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  ggtitle(paste0("Percent Positive HHV8 Swabs for ",desc," Study ARM")) +
  scale_color_manual(values=colors, name="Shedding %",aesthetics = c("fill")) +
  ylim(0,100) +
  xlab("PTIDs") + ylab("Percent Positive")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

dev.off()
