library(RColorBrewer)
# Load ggplot2
library(ggplot2)
library(tidyverse)

args<-commandArgs(trailingOnly=T)
ptid_list<-args[1]

ptids  <- file(ptid_list, open = "r")

names=list();
sim1_shed_perc=list();
sim2_shed_perc=list();
sim1_med_hhv8=list();
sim2_med_hhv8=list();
first_line=1;

while (length(oneLine <- readLines(ptids, n = 1, warn = FALSE)) > 0) {

    ptid = oneLine;

    sim_data = read.csv(paste0(ptid,".dat1.csv"))
    #time,vet1,vet2,vet3,vet4,vet5,vet6,vet7,vet8,vet9,vet10

    sim1_hhv8 = 0
    for (i in 1:length(sim_data$vet1)) {
	if (sim_data$vet1[i] > 0) {
	    sim1_hhv8 = sim1_hhv8+1;
	}
    }
    sim1_hhv8_vl=median(sim_data$vet1)
    sim1_hhv8_shed=100*sim1_hhv8 / length(sim_data$vet1)

    sim2_hhv8 = 0
    for (i in 1:length(sim_data$vet2)) {
	if (sim_data$vet2[i] > 0) {
	    sim2_hhv8 = sim2_hhv8+1;
	}
    }
    sim2_hhv8_vl=median(sim_data$vet2)
    sim2_hhv8_shed=100*sim2_hhv8 / length(sim_data$vet2)

    ptid_num = sub("ptid_", "", ptid)
    if (first_line==1) {
	sim1_shed_perc = sim1_hhv8_shed
	sim2_shed_perc = sim2_hhv8_shed
	sim1_med_hhv8 = sim1_hhv8_vl
	sim2_med_hhv8 = sim2_hhv8_vl
	print(paste("first",ptid))
    } else {
	sim1_shed_perc= c(sim1_shed_perc,sim1_hhv8_shed)
	sim2_shed_perc= c(sim2_shed_perc,sim2_hhv8_shed)
	sim1_med_hhv8= c(sim1_med_hhv8,sim1_hhv8_vl)
	sim2_med_hhv8= c(sim2_med_hhv8,sim2_hhv8_vl)
	print(paste("added",ptid))
    }
    first_line=0
}

pdf("ptid_shed_perc3.pdf",height=8,width=8)

cor.test(sim2_shed_perc, sim1_shed_perc)

df<-data.frame(x=sim2_shed_perc, y=sim1_shed_perc)

ggplot(df, aes(x=x, y=y)) +
  geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Shedding Perc Sim Run 1 vs 2") +
  xlab("Run 1") + ylab("Run 2") +
  xlim(0, 100) +
  ylim(0, 100)

dev.off()

pdf("ptid_med_hhv8_vl3.pdf",height=8,width=8)

cor.test(sim2_med_hhv8, sim1_med_hhv8)

df<-data.frame(x=sim2_med_hhv8, y=sim1_med_hhv8)

ggplot(df, aes(x=x, y=y)) +
  geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Med Log HHV8 Sim Run 1 vs 2") +
  xlab("Run 1") + ylab("Run 2") +
  xlim(0, 8) +
  ylim(0, 8)

dev.off()
