library(tidyverse)
library(viridis)
library(dplyr)
library(svglite)


SY_qPCR <- read.csv("SY5Y_Tidy_4.csv")
SY_qPCR$Treatment <- factor(SY_qPCR$Treatment,levels=c("DMSO","CHX"))

agg=aggregate(M_Value~Treatment*Gene, data=SY_qPCR, FUN="mean")
standard_error <- function(x) sd(x) / sqrt(length(x))
agg$se=aggregate(M_Value~Treatment*Gene, data=SY_qPCR, FUN="standard_error")$M_Value

limits <- aes(ymax=M_Value+se, ymin=M_Value-se)
dodge <- position_dodge(width=0.9) 

p <- ggplot(agg, aes(fill=Treatment, y=M_Value, x=Gene)) 

p+geom_bar(position=dodge, stat="identity") +
  xlab("") +
  ylab("") +
  ggtitle("SH-SY5Y_qPCR") +
  geom_errorbar(limits, position=dodge, width=0.4)+
  geom_point(data=SY_qPCR, aes(x=Gene), shape=20, position = 
               position_jitterdodge(jitter.width = 0.5, jitter.height=0, 
                                    dodge.width=0.9))+
  theme_bw()+
  theme(title=element_text(face="bold"))+
  theme(legend.title=element_text(size=14)) +
  theme(legend.text=element_text(size=12,face="bold")) +
  theme(axis.text=element_text(color="black",size=12,face="bold")) +
  theme(axis.title.y=element_text(color="black",size=18,face="bold")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  scale_fill_viridis_d(begin=0.3, end=1)




ggsave(file="SY5Y_qPCR_Matt_Paper.pdf", device = "pdf", width=5, height=4)
ggsave(file="SY5Y_qPCR_Matt_Paper.ps", device = "ps", width=5, height=4)
ggsave(file="SY5Y_qPCR_Matt_Paper.svg", device = "svg", width=5, height=4)

ps_SY.HNRNPL <- shapiro.test(filter(SY_qPCR, Treatment == "CHX" & Gene == "HNRNPL")$Value)
ps_SY.UNC13A <- shapiro.test(filter(SY_qPCR, Treatment == "CHX" & Gene == "UNC13A")$Value)
ps_SY.UNC13B <- shapiro.test(filter(SY_qPCR, Treatment == "CHX" & Gene == "UNC13B")$M_Value)
ps_SY.STMN2 <- shapiro.test(filter(SY_qPCR, Treatment == "CHX" & Gene == "STMN2")$Value)
#Data passes normality test
p_SY.HNRNPL <- t.test(filter(SY_qPCR, Treatment == "CHX" & Gene == "HNRNPL")$M_Value, mu=100)
p_SY.UNC13A <- t.test(filter(SY_qPCR, Treatment == "CHX" & Gene == "UNC13A")$M_Value, mu=100)
p_SY.UNC13B <- t.test(filter(SY_qPCR, Treatment == "CHX" & Gene == "UNC13B")$M_Value, mu=100)
p_SY.UNC13B_Avg <- t.test(filter(SY_qPCR, Treatment == "CHX" & Gene == "UNC13B")$Avg, mu=100)
p_SY.STMN2 <- t.test(filter(SY_qPCR, Treatment == "CHX" & Gene == "STMN2")$M_Value, mu=100)

p_SY.TDP <- t.test(filter(SY_qPCR, Treatment == "shTDP-43" & Gene == "TDP-43")$Value, mu=100)
p_SY.13A <- t.test(filter(SY_qPCR, Treatment == "shTDP-43" & Gene == "UNC13A")$Value, mu=100)
p_SY.13B <- t.test(filter(SY_qPCR, Treatment == "shTDP-43" & Gene == "UNC13B")$Value, mu=100)

# TDP-43 p = 1.3 e-5, UNC13A p = 7.6 e-5. UNC13B p = 0.002

