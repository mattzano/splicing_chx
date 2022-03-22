# NMD experiment on SH-SY5Y cells - CHX
library(rstatix)
library(ggnewscale)
library(patchwork)
std <- function(x) sd(x)/sqrt(length(x))

for_bar_plot_input <- fread(file.path(here::here(),"data","SY5Y_Tidy_4.csv"))
for_bar_plot <- data.frame(for_bar_plot_input)
for_bar_plot$Treatment <- factor(for_bar_plot$Treatment, levels = c("DMSO", "CHX"))
shapiro.test(filter(for_bar_plot, Treatment == "CHX" & Gene == "HNRNPL")$Value)
shapiro.test(filter(for_bar_plot, Treatment == "CHX" & Gene == "STMN2")$Value)
shapiro.test(filter(for_bar_plot, Treatment == "CHX" & Gene == "UNC13A")$Value)
shapiro.test(filter(for_bar_plot, Treatment == "CHX" & Gene == "UNC13B")$Value)
stattest <- for_bar_plot %>%
  group_by(Gene) %>%
  pairwise_t_test(Value ~ Treatment, ref.group = "DMSO", alternative = "greater") %>%
  add_significance() %>% 
  add_xy_position(x="Treatment")
stattest
sum_for_bar_plot <- for_bar_plot %>% 
  group_by(Gene,Treatment) %>% 
  summarize(mean = mean(Value),sem = std(Value)) 
ggplot(sum_for_bar_plot) +
  geom_col(aes(x = Gene, fill = Treatment, y = mean),position = "dodge2",color = 'black') +
  scale_fill_manual(values = c("#F9AC66","#B4C640")) +
  new_scale_fill() +
  geom_errorbar(aes(x = Gene, ymin = mean - sem, ymax = mean + sem),
                position = "dodge2",size = 1.8) +
  geom_point(data = for_bar_plot, aes(x = Gene, fill = Treatment, y = Value),
             pch = 21, 
             size = 4, stroke = 1.5,
             position = position_jitterdodge(jitter.width = 0.2)) +
  scale_fill_manual(values = c("#000000","#000000")) +
  ggpubr::theme_pubr() + 
  theme(axis.ticks.length=unit(0.1,"inch"),
        axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 2))

# NMD experiment on SH-SY5Y cells - sTDP43+sUPF1

for_bar_plot_sy_input <- fread(file.path(here::here(),"data","results_upf1.csv"))
for_bar_plot_sy_df <- data.frame(for_bar_plot_sy_input)
for_bar_plot_sy_df$Condition <- as.factor(for_bar_plot_sy_df$Condition)
for_bar_plot_sy_df$Value <- 100 * for_bar_plot_sy_df$Value
#assess TDP-43 and UPF1 KD
for_bar_plot_sz <- for_bar_plot_sy_df %>%
  filter(Gene == "TDP43" & (Condition == "siTDP43+siCTRL" | Condition == "siTDP43+siUPF1" | Condition == "Non-treated"))
shapiro.test(filter(for_bar_plot_sz, Condition == "siTDP43+siCTRL" & Gene == "TDP43")$Value)
shapiro.test(filter(for_bar_plot_sz, Condition == "siTDP43+siUPF1" & Gene == "TDP43")$Value)
stattest_sz <- for_bar_plot_sz %>%
  pairwise_t_test(Value ~ Condition, ref.group = "Non-treated", alternative = "less") %>%
  add_significance() %>% 
  add_xy_position(x="Condition")
stattest_sz
for_bar_plot_sx <- for_bar_plot_sy_df %>%
  filter(Gene == "UPF1" & (Condition == "siTDP43+siCTRL" | Condition == "siTDP43+siUPF1"))
shapiro.test(filter(for_bar_plot_sx, Condition == "siTDP43+siUPF1" & Gene == "UPF1")$Value)
stattest_sx <- for_bar_plot_sx %>%
  pairwise_t_test(Value ~ Condition, ref.group = "siTDP43+siCTRL", alternative = "less") %>%
  add_significance() %>% 
  add_xy_position(x="Condition")
stattest_sx
for_bar_plot_szx <- for_bar_plot_sy_df %>%
  filter(Gene == "TDP43" & (Condition == "siTDP43+siCTRL" | Condition == "siTDP43+siUPF1" | Condition == "Non-treated")
         | Gene == "UPF1" & (Condition == "siTDP43+siCTRL" | Condition == "siTDP43+siUPF1"))
#fix position of stats for plot
stattest_sx$xmin[1] <- 1.8
stattest_sx$xmax[1] <- 2.2
stattest_sx$y.position[1] <- 1.3
stattest_sz$xmin[1] <- 0.75
stattest_sz$xmin[2] <- 0.75
stattest_sz$xmax[1] <- 1
stattest_sz$xmax[2] <- 1.25
stattest_sz$y.position[1] <- 1.2
stattest_sz$y.position[2] <- 1.3
nmd_plotz = ggbarplot(for_bar_plot_szx,
                      x = 'Gene',
                      add = c("mean_se","jitter"),
                      y = 'Value',
                      color = 'Condition',
                      position = position_dodge(0.8),
                      dot.size = 10) +
  scale_y_continuous() + ylab("Percent of transcript after UPF1 knockdown") + xlab("Gene") + 
  ggpubr::theme_pubr() + stat_pvalue_manual(stattest_sz) + stat_pvalue_manual(stattest_sx)
plot(nmd_plotz)
#assess rescue of NMD-sensitive transcripts after UPF1 inhibition
for_bar_plot_sy <- for_bar_plot_sy_df %>%
  filter(Gene == "hnRNPL" | Gene == "STMN2" | Gene == "UNC13A" | Gene == "UNC13B")
shapiro.test(filter(for_bar_plot_sy, Condition == "siTDP43+siUPF1" & Gene == "hnRNPL")$Value)
shapiro.test(filter(for_bar_plot_sy, Condition == "siTDP43+siUPF1" & Gene == "STMN2")$Value)
shapiro.test(filter(for_bar_plot_sy, Condition == "siTDP43+siUPF1" & Gene == "UNC13A")$Value)
shapiro.test(filter(for_bar_plot_sy, Condition == "siTDP43+siUPF1" & Gene == "UNC13B")$Value)
stattest_sy <- for_bar_plot_sy %>%
  group_by(Gene) %>%
  pairwise_t_test(Value ~ Condition, ref.group = "siTDP43+siCTRL", alternative = "greater") %>%
  add_significance() %>% 
  add_xy_position(x="Condition")
stattest_sy

nmd_plot = ggbarplot(for_bar_plot_sy, 
                     x = 'Gene', 
                     add = c("mean_se","jitter"), 
                     y = 'Value',
                     color = 'Condition',
                     fill = "Condition",
                     position = position_dodge(0.8), 
                     dot.size = 20,
                     size = 1.5) +
  scale_y_continuous() + 
  ggpubr::theme_pubr() + ylab("Percent of transcript after UPF1 knockdown") + xlab("Gene") +
  stat_pvalue_manual(stattest_sy, label="p.signif", x = "Gene") +
  theme(axis.ticks.length=unit(0.1,"inch"),
        axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 2))+
  scale_fill_manual(values = c("#fca361ff","#aac043ff")) +
  scale_color_manual(values = c("#000000ff","#000000ff"))
plot(nmd_plot)

sum_for_bar_plot_sy <- for_bar_plot_sy %>% 
  group_by(Gene,Condition) %>% 
  summarize(mean = mean(Value),sem = std(Value)) 

sy_plot <- ggplot(sum_for_bar_plot_sy) +
  geom_col(aes(x = Gene, y = mean,
          #add = c("mean_se","jitter"), 
          fill = Condition),
          position = "dodge2", 
          #dot.size = 20,
          #size = 1.5,
          color = 'black') +
  scale_fill_manual(values = c("#F9AC66","#B4C640")) +
  scale_y_continuous() + 
  new_scale_fill() +
  geom_errorbar(aes(x = Gene, ymin = mean - sem, ymax = mean + sem),
                position = "dodge2", size = 1.1) +
  geom_point(data = for_bar_plot_sy, aes(x = Gene, fill = Condition, y =Value),
             pch = 21, 
             size = 2, stroke = 0.8,
             position = position_jitterdodge(jitter.width = 0.2)) +
  scale_fill_manual(values = c("#000000","#000000")) +
  ggpubr::theme_pubr() + 
  ylab("Percent of transcript after UPF1 KD") + xlab("Gene") +
  
  stat_pvalue_manual(stattest_sy, label="p.signif", x = "Gene", size = 6) +
  theme(axis.ticks.length=unit(0.1,"inch"),
        axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 2))
  
plot <- ggplot(sum_for_bar_plot) +
  geom_col(aes(x = Gene, fill = Treatment, y = mean),position = "dodge2",color = 'black') +
  scale_fill_manual(values = c("#F9AC66","#B4C640")) +
  new_scale_fill() +
  geom_errorbar(aes(x = Gene, ymin = mean - sem, ymax = mean + sem),
                position = "dodge2",size = 1.1) +
  geom_point(data = for_bar_plot, aes(x = Gene, fill = Treatment, y = Value),
             pch = 21, 
             size = 2, stroke = 0.8,
             position = position_jitterdodge(jitter.width = 0.2)) +
  scale_fill_manual(values = c("#000000","#000000")) +
  ggpubr::theme_pubr() + 
  stat_pvalue_manual(stattest, label="p.signif", x = "Gene", size = 6) +
  ylab("Percent of transcript after CHX treament") + xlab("") +
  scale_x_discrete(labels = c("", "", "", "")) +
  theme(axis.ticks.length=unit(0.1,"inch"),
        axis.line = element_line(colour = 'black', size = 1.5),
        axis.ticks = element_line(colour = "black", size = 2))

sum_plot <- plot / sy_plot
plot(sum_plot)
ggsave("patch_nmd.png", sum_plot)
