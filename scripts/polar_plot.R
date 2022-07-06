####polar_plot
data <- as.data.frame(table(big_changing_disease_labelled2$Level3,
                            big_changing_disease_labelled2$category))

colnames(data) <- c("individual", "group", "value")
data$relfreqz  <- ave(data$value, data$group, FUN=function(x) x/sum(x))

big_disease <- rbindlist(datalist_disease, idcol = TRUE)
big_disease$.id <- factor(big_disease$.id, levels = input_disease)
big_disease_table <- as.data.frame(table(big_disease$Level3))
colnames(big_disease_table) <- c("individual", "value")
data <- data %>% left_join(big_disease_table, by = "individual")

data$relfreq  <- (data$value.x/data$value.y) *1000


empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=relfreq, fill=individual)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=relfreq, fill=individual), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 150, xend = start, yend = 150), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(50, 100, 150, 200), label = c("5%", "10%", "15%", "20%") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=relfreq, fill=individual), stat="identity", alpha=0.5) +
  ylim(-100,300) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=relfreq+10, label=individual, hjust=hjust), colour="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -40, label=group), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

p





data <- data %>%
  #group_by(individual) %>%
  #dplyr::filter(value > 10) #%>%
  dplyr::mutate(new_group = ifelse(group == "BIALLELIC" | group == "BOTH" | group == "X-LINKED", "Loss of Function",
                                   ifelse(group == "MONOALLELIC", "Gain of function", "Other")))

data$relfreqz  <- ave(data$value, data$individual, FUN=function(x) x/sum(x))

  #group_by(individual) %>%
  #dplyr::filter(n() > 10) %>%
data %>%
  ggplot(aes(x = factor(individual), y = relfreqz, fill = new_group)) +
  geom_col(position=position_stack(), width = .7) +
  #geom_text(aes(label = round(relfreqz, 2)), position = position_stack(vjust = 0.5), size = 2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

cut_prop = data %>% 
  group_by(individual) %>% 
  summarise(cut_prop = n()/nrow(data))

diamonds = left_join(data, cut_prop)

df <- diamonds %>% 
  dplyr::count(individual, new_group, relfreqz, cut_prop) %>% 
  group_by(individual) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup

ggplot(data = df,
       aes(x = individual, fill = new_group, y = relfreqz, width = cut_prop)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip()


big_changing_mono <- big_changing_disease_labelled %>%
  dplyr::filter(category == "MONOALLELIC") %>%
  dplyr::filter(junc_type != "annotated" & junc_type != "none") %>%
  distinct(paste_into_igv_junction, .keep_all = T)


big_delta_mono <- big_delta_disease_labelled %>%
  dplyr::filter(category == "MONOALLELIC") %>%
  dplyr::filter(junc_cat != "annotated" & junc_cat != "none") %>%
  distinct(paste_into_igv_junction, .keep_all = T)


big_cryptic_mono <- big_cryptic_disease_labelled %>%
  dplyr::filter(category == "MONOALLELIC") %>%
  dplyr::filter(junc_type != "annotated" & junc_type != "none") %>%
  distinct(paste_into_igv_junction, .keep_all = T)
