library(RColorBrewer)
# Load ggplot2
library(ggplot2)
library(tidyverse)

args<-commandArgs(trailingOnly=T)
ptid_list<-args[1]
desc<-args[2]

ptids  <- file(ptid_list, open = "r")

names=list();
cols=list();
vets=list();

while (length(oneLine <- readLines(ptids, n = 1, warn = FALSE)) > 0) {

    ptid = oneLine;

    ptid_data = read.csv(paste0(ptid,".dat1.csv"))
    #time,vet1,vet2,vet3,vet4,vet5,vet6,vet7,vet8,vet9,vet10

    hhv8 = 0
    for (i in 1:length(ptid_data$vet1)) {
	if (ptid_data$vet1[i] > 0) {
	    hhv8 = hhv8+1;
	}
    }
    #hhv8_stat=median(ptid_data$vet1)
    hhv8_stat=100*hhv8 / length(ptid_data$vet1)
    print(hhv8_stat)

    if (100*hhv8 / length(ptid_data$vet1) < 20) {
	new_col="<20"
    } else if (100*hhv8 / length(ptid_data$vet1) < 40) {
	new_col="<40"
    } else if (100*hhv8 / length(ptid_data$vet1) < 60) {
	new_col="<60"
    } else if (100*hhv8 / length(ptid_data$vet1) < 80) {
	new_col="<80"
    } else {
	new_col=">=80"
    }
    ptid_num = sub("ptid_", "", ptid)
    if (length(cols)==0) {
	vl_stats = rep(hhv8_stat,length(ptid_data$vet1))
	names = rep(ptid_num,length(ptid_data$vet1))
	cols = rep(new_col,length(ptid_data$vet1))
	vets = ptid_data$vet1
	print(paste("first",ptid))
    } else {
	vl_stats= c(vl_stats,rep(hhv8_stat,length(ptid_data$vet1)))
	names= c(names,rep(ptid_num,length(ptid_data$vet1)))
	cols= c(cols,rep(new_col,length(ptid_data$vet1)))
	vets = c(vets,ptid_data$vet1)
	print(paste("added",ptid))
    }
}

y_lim=c(0,8)

pdf(paste0("ptid_",desc,"_shedders.pdf"),height=8,width=16)

df<-data.frame(vals=vets, names=names, fill=cols, vl_stat=vl_stats)

sorted_df <- df[order(df$vl_stat),]
options(max.print=999999)
print(sorted_df)
  
colors <- c("<20" = "lightblue", "<40" = "blue", "<60" = "yellow", "<80" = "orange", ">=80" = "red")

ggplot(sorted_df, aes(x=fct_inorder(names), y=vals, fill=fill)) +
  geom_boxplot() +
  ggtitle(paste0("PTID HHV8 ",desc," Shedders")) +
  scale_color_manual(values=colors, name="Shedding %",aesthetics = c("fill")) +
  ylim(0,8) +
  xlab("PTID") + ylab("log HHV8")

dev.off()
