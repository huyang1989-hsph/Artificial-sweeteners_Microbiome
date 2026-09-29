library(meta)
library(ggplotify)
library(patchwork)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyverse)
###########################################################################################################################
df.all <- read.csv("~/replicate_bb4_score_or.csv",header=T)

df.all$study<-ifelse(df.all$study=="pooled","Pooled",df.all$study)

# Create flag for truncated CIs
df.all <- df.all %>%
  mutate(
    study = factor(study, levels = rev(unique(study))),
    score = factor(score, levels = unique(score)),
    weight = 1 / se,
    is_pooled = ifelse(study == "Pooled", "Pooled", "Study"),
    shape_group = ifelse(is_pooled == "Pooled", "Diamond", "Square"),
    uci_capped = ifelse(uci > 2.5, 2.5, uci),
    lci_capped = ifelse(lci < 0.4, 0.4, lci),
    over_limit = uci > 2.5,
    under_limit = lci < 0.4,
    pooled_label = ifelse(
      study == "Pooled",
      sprintf("%.2f(%.2f,%.2f)", or, lci, uci),
      NA
    ),
    sex="All"
    )

df.all.ss<-df.all %>% filter(score %in% c("Acesulfame K", 
                                  "Aspartame", 
                                  "Saccharin", 
                                  "Sucralose"))
df.all.sg<-df.all %>% filter(score %in% c("Erythritol",
                                  "Maltitol",
                                  "Mannitol",
                                  "Sorbitol",
                                  "Xylitol"))
df.all.ssg<-df.all %>% filter(score %in% c("Sugar alcohols",
                                   "Synthetic sweeteners"))

fig3.1<-
  ggplot(df.all.ss, aes(x = or, y = study, color = is_pooled)) +
  geom_point(aes(size = weight, shape = shape_group), alpha = 0.9)  +
  geom_errorbarh(aes(xmin = lci_capped, xmax = uci_capped), height = 0.1) +
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
 
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Study" = "black", "Pooled" = "red"), guide = "none") +
  # coord_cartesian(xlim = c(0, 2.5)) +  # Set lower and upper x limits
  facet_wrap(~score, scales = "free_y", ncol = 4) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  
  geom_text(
    data = df.all.ss %>% filter(study == "Pooled"),
    aes(label = pooled_label),
    hjust = 0.5,  # push label slightly right of point
    vjust= 2,
    size = 6,
    color = "red"
  )+
  theme_minimal() +
  scale_x_log10() +
  theme(

    legend.position = "none",               # Remove legend
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_text(size = 12, color = "black",face = "bold"),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    axis.title.x = element_text(size = 15, color = "black",face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing.x = unit(0.5, "cm")
  )

png("~/fig3_1_all.png",res=350,width=9000, height=4000)
fig3.1
dev.off()


fig3.2<-
  ggplot(df.all.sg, aes(x = or, y = study, color = is_pooled)) +
  geom_point(aes(size = weight, shape = shape_group), alpha = 0.9)  +
  geom_errorbarh(aes(xmin = lci_capped, xmax = uci_capped), height = 0.1) +
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
 
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Study" = "black", "Pooled" = "red"), guide = "none") +
  # coord_cartesian(xlim = c(0, 2.5)) +  # Set lower and upper x limits
  facet_wrap(~score, scales = "free_y", ncol = 5) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  
  geom_text(
    data = df.all.sg %>% filter(study == "Pooled"),
    aes(label = pooled_label),
    hjust = 0.5,  # push label slightly right of point
    vjust= 2,
    size = 6,
    color = "red"
  )+
  theme_minimal() +
  scale_x_log10() +
  theme(
    
    legend.position = "none",               # Remove legend
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_text(size = 12, color = "black",face = "bold"),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    axis.title.x = element_text(size = 15, color = "black",face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing.x = unit(0.5, "cm")
  )

png("~/fig3_2_all.png",res=350,width=12000, height=4000)
fig3.2
dev.off()


fig3.3<-
  ggplot(df.all.ssg, aes(x = or, y = study, color = is_pooled)) +
  geom_point(aes(size = weight, shape = shape_group), alpha = 0.9)  +
  geom_errorbarh(aes(xmin = lci_capped, xmax = uci_capped), height = 0.1) +
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Study" = "black", "Pooled" = "red"), guide = "none") +
  # coord_cartesian(xlim = c(0, 2.5)) +  # Set lower and upper x limits
  facet_wrap(~score, scales = "free_y", ncol = 2) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  
  geom_text(
    data = df.all.ssg %>% filter(study == "Pooled"),
    aes(label = pooled_label),
    hjust = 0.5,  # push label slightly right of point
    vjust= 2,
    size = 6,
    color = "red"
  )+
  theme_minimal() +
  scale_x_log10() +
  theme(
    
    legend.position = "none",               # Remove legend
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_text(size = 12, color = "black",face = "bold"),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    axis.title.x = element_text(size = 15, color = "black",face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

png("~/fig3_3_all.png",res=350,width=8000, height=4000)
fig3.3
dev.off()
###########################################################################################################################
#################################### MALE #################################### 
###########################################################################################################################
df.men <- read.csv("~/replicate_bb4_score_male.csv",header=T)

df.men$study<-ifelse(df.men$study=="pooled","Pooled",df.men$study)

# Create flag for truncated CIs
df.men <- df.men %>%
  mutate(
    study = factor(study, levels = rev(unique(study))),
    score = factor(score, levels = unique(score)),
    weight = 1 / se,
    is_pooled = ifelse(study == "Pooled", "Pooled", "Study"),
    shape_group = ifelse(is_pooled == "Pooled", "Diamond", "Square"),
    uci_capped = ifelse(uci > 3, 3, uci),
    lci_capped = ifelse(lci < 0.2, 0.2, lci),
    over_limit = uci > 3,
    under_limit = lci < 0.2,
    pooled_label = ifelse(
      study == "Pooled",
      sprintf("%.2f(%.2f,%.2f)", or, lci, uci),
      NA
    ),
    sex="Men"
  )

df.ss<-df %>% filter(score %in% c("Acesulfame K", 
                                  "Aspartame", 
                                  "Saccharin", 
                                  "Sucralose"))
df.sg<-df %>% filter(score %in% c("Erythritol",
                                  "Maltitol",
                                  "Mannitol",
                                  "Sorbitol",
                                  "Xylitol"))
df.ssg<-df %>% filter(score %in% c("Sugar alcohols",
                                   "Synthetic sweeteners"))

fig3.1<-
  ggplot(df.ss, aes(x = or, y = study, color = is_pooled)) +
  geom_point(aes(size = weight, shape = shape_group), alpha = 0.9)  +
  geom_errorbarh(aes(xmin = lci_capped, xmax = uci_capped), height = 0.1) +
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
  # Right-pointing arrows for uci > 2
  geom_segment(
    data = df.ss %>% filter(over_limit),
    aes(x = uci_capped, xend = 3.05, y = study, yend = study),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    color = "black"
  ) +

  # Left-pointing arrows for lci < 0.6
  geom_segment(
    data = df.ss %>% filter(under_limit),
    aes(x = lci_capped, xend = 0.15, y = study, yend = study),  # x > xend
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),  # default ends = "last"
    color = "black"
  ) +

  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Study" = "black", "Pooled" = "red"), guide = "none") +
  
  # Customize x-axis
  scale_x_continuous(
    name = "Odds Ratio",  # Change axis title
    breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)
  ) +
  
  coord_cartesian(xlim = c(0, 3)) +  # Set lower and upper x limits
  facet_wrap(~score, scales = "free_y", ncol = 4) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  
  geom_text(
    data = df.ss %>% filter(study == "Pooled"),
    aes(label = pooled_label),
    hjust = 0.5,  # push label slightly right of point
    vjust= 2,
    size = 6,
    color = "red"
  )+
  theme_minimal() +
  theme(
    
    legend.position = "none",               # Remove legend
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    axis.title.x = element_text(size = 15, color = "black",face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing.x = unit(0.5, "cm")
  )

png("~/fig3_1_male.png",res=350,width=9000, height=4000)
fig3.1
dev.off()


fig3.2<-
  ggplot(df.sg, aes(x = or, y = study, color = is_pooled)) +
  geom_point(aes(size = weight, shape = shape_group), alpha = 0.9)  +
  geom_errorbarh(aes(xmin = lci_capped, xmax = uci_capped), height = 0.1) +
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
  # Right-pointing arrows for uci > 2
  geom_segment(
    data = df.sg %>% filter(over_limit),
    aes(x = uci_capped, xend = 3.05, y = study, yend = study),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    color = "black"
  ) +

  # Left-pointing arrows for lci < 0.6
  geom_segment(
    data = df.sg %>% filter(under_limit),
    aes(x = lci_capped, xend = 0.15, y = study, yend = study),  # x > xend
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),  # default ends = "last"
    color = "black"
  ) +

  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Study" = "black", "Pooled" = "red"), guide = "none") +
  
  # Customize x-axis
  scale_x_continuous(
    name = "Odds Ratio",  # Change axis title
    breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)
  ) +
  
  coord_cartesian(xlim = c(0, 3)) +  # Set lower and upper x limits
  facet_wrap(~score, scales = "free_y", ncol = 5) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  
  geom_text(
    data = df.sg %>% filter(study == "Pooled"),
    aes(label = pooled_label),
    hjust = 0.5,  # push label slightly right of point
    vjust= 2,
    size = 6,
    color = "red"
  )+
  theme_minimal() +
  theme(
    
    legend.position = "none",               # Remove legend
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    axis.title.x = element_text(size = 15, color = "black",face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing.x = unit(0.5, "cm")
  )

png("~/fig3_2_male.png",res=350,width=12000, height=4000)
fig3.2
dev.off()


fig3.3<-
  ggplot(df.ssg, aes(x = or, y = study, color = is_pooled)) +
  geom_point(aes(size = weight, shape = shape_group), alpha = 0.9)  +
  geom_errorbarh(aes(xmin = lci_capped, xmax = uci_capped), height = 0.1) +
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled

  # Right-pointing arrows for uci > 2
  geom_segment(
    data = df.ssg %>% filter(over_limit),
    aes(x = uci_capped, xend = 3.03, y = study, yend = study),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    color = "black"
  ) +
  # Left-pointing arrows for lci < 0.6
  geom_segment(
    data = df.ssg %>% filter(under_limit),
    aes(x = lci_capped, xend = 0.28, y = study, yend = study),  # x > xend
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),  # default ends = "last"
    color = "black"
  ) +
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Study" = "black", "Pooled" = "red"), guide = "none") +
  
  # Customize x-axis
  scale_x_continuous(
    name = "Odds Ratio",  # Change axis title
    breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)
  ) +
  
  coord_cartesian(xlim = c(0, 3)) +  # Set lower and upper x limits
  facet_wrap(~score, scales = "free_y", ncol = 2) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  
  geom_text(
    data = df.ssg %>% filter(study == "Pooled"),
    aes(label = pooled_label),
    hjust = 0.5,  # push label slightly right of point
    vjust= 2,
    size = 6,
    color = "red"
  )+
  theme_minimal() +
  theme(
    
    legend.position = "none",               # Remove legend
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    axis.title.x = element_text(size = 15, color = "black",face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing.x = unit(0.5, "cm")
  )

png("~/fig3_3_male.png",res=350,width=8000, height=4000)
fig3.3
dev.off()
###########################################################################################################################
#################################### MALE #################################### 
###########################################################################################################################
df.women <- read.csv("~/replicate_bb4_score_female.csv",header=T)

df.women$study<-ifelse(df.women$study=="pooled_no_directplus","Pooled",df.women$study)

# Create flag for truncated CIs
df.women <- df.women %>%
  mutate(
    study = factor(study, levels = rev(unique(study))),
    score = factor(score, levels = unique(score)),
    weight = 1 / se,
    is_pooled = ifelse(study == "Pooled", "Pooled", "Study"),
    shape_group = ifelse(is_pooled == "Pooled", "Diamond", "Square"),
    uci_capped = ifelse(uci > 3, 3, uci),
    lci_capped = ifelse(lci < 0.2, 0.2, lci),
    over_limit = uci > 3,
    under_limit = lci < 0.2,
    pooled_label = ifelse(
      study == "Pooled",
      sprintf("%.2f(%.2f,%.2f)", or, lci, uci),
      NA
    ),
  sex="Women"
  )

df.ss<-df %>% filter(score %in% c("Acesulfame K", 
                                  "Aspartame", 
                                  "Saccharin", 
                                  "Sucralose"))
df.sg<-df %>% filter(score %in% c("Erythritol",
                                  "Maltitol",
                                  "Mannitol",
                                  "Sorbitol",
                                  "Xylitol"))
df.ssg<-df %>% filter(score %in% c("Sugar alcohols",
                                   "Synthetic sweeteners"))

fig3.1<-
  ggplot(df.ss, aes(x = or, y = study, color = is_pooled)) +
  geom_point(aes(size = weight, shape = shape_group), alpha = 0.9)  +
  geom_errorbarh(aes(xmin = lci_capped, xmax = uci_capped), height = 0.1) +
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
  # Right-pointing arrows for uci > 2
  geom_segment(
    data = df.ss %>% filter(over_limit),
    aes(x = uci_capped, xend = 3.05, y = study, yend = study),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    color = "black"
  ) +
  
  # Left-pointing arrows for lci < 0.6
  geom_segment(
    data = df.ss %>% filter(under_limit),
    aes(x = lci_capped, xend = 0.15, y = study, yend = study),  # x > xend
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),  # default ends = "last"
    color = "black"
  ) +
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Study" = "black", "Pooled" = "red"), guide = "none") +
  
  # Customize x-axis
  scale_x_continuous(
    name = "Odds Ratio",  # Change axis title
    breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)
  ) +
  
  coord_cartesian(xlim = c(0, 3)) +  # Set lower and upper x limits
  facet_wrap(~score, scales = "free_y", ncol = 4) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  
  geom_text(
    data = df.ss %>% filter(study == "Pooled"),
    aes(label = pooled_label),
    hjust = 0.5,  # push label slightly right of point
    vjust= 2,
    size = 6,
    color = "red"
  )+
  theme_minimal() +
  theme(
    
    legend.position = "none",               # Remove legend
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    axis.title.x = element_text(size = 15, color = "black",face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing.x = unit(0.5, "cm")
  )

png("~/fig3_1_female.png",res=350,width=9000, height=4000)
fig3.1
dev.off()


fig3.2<-
  ggplot(df.sg, aes(x = or, y = study, color = is_pooled)) +
  geom_point(aes(size = weight, shape = shape_group), alpha = 0.9)  +
  geom_errorbarh(aes(xmin = lci_capped, xmax = uci_capped), height = 0.1) +
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
  # Right-pointing arrows for uci > 2
  geom_segment(
    data = df.sg %>% filter(over_limit),
    aes(x = uci_capped, xend = 3.05, y = study, yend = study),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    color = "black"
  ) +
  
  # Left-pointing arrows for lci < 0.6
  geom_segment(
    data = df.sg %>% filter(under_limit),
    aes(x = lci_capped, xend = 0.15, y = study, yend = study),  # x > xend
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),  # default ends = "last"
    color = "black"
  ) +
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Study" = "black", "Pooled" = "red"), guide = "none") +
  
  # Customize x-axis
  scale_x_continuous(
    name = "Odds Ratio",  # Change axis title
    breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)
  ) +
  
  coord_cartesian(xlim = c(0, 3)) +  # Set lower and upper x limits
  facet_wrap(~score, scales = "free_y", ncol = 5) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  
  geom_text(
    data = df.sg %>% filter(study == "Pooled"),
    aes(label = pooled_label),
    hjust = 0.5,  # push label slightly right of point
    vjust= 2,
    size = 6,
    color = "red"
  )+
  theme_minimal() +
  theme(
    
    legend.position = "none",               # Remove legend
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    axis.title.x = element_text(size = 15, color = "black",face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing.x = unit(0.5, "cm")
  )

png("~/fig3_2_female.png",res=350,width=12000, height=4000)
fig3.2
dev.off()


fig3.3<-
  ggplot(df.ssg, aes(x = or, y = study, color = is_pooled)) +
  geom_point(aes(size = weight, shape = shape_group), alpha = 0.9)  +
  geom_errorbarh(aes(xmin = lci_capped, xmax = uci_capped), height = 0.1) +
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
  
  # Right-pointing arrows for uci > 2
  geom_segment(
    data = df.ssg %>% filter(over_limit),
    aes(x = uci_capped, xend = 3.03, y = study, yend = study),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    color = "black"
  ) +
  # Left-pointing arrows for lci < 0.6
  geom_segment(
    data = df.ssg %>% filter(under_limit),
    aes(x = lci_capped, xend = 0.28, y = study, yend = study),  # x > xend
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),  # default ends = "last"
    color = "black"
  ) +
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Study" = "black", "Pooled" = "red"), guide = "none") +
  
  # Customize x-axis
  scale_x_continuous(
    name = "Odds Ratio",  # Change axis title
    breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)
  ) +
  
  coord_cartesian(xlim = c(0, 3)) +  # Set lower and upper x limits
  facet_wrap(~score, scales = "free_y", ncol = 2) +
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  
  geom_text(
    data = df.ssg %>% filter(study == "Pooled"),
    aes(label = pooled_label),
    hjust = 0.5,  # push label slightly right of point
    vjust= 2,
    size = 6,
    color = "red"
  )+
  theme_minimal() +
  theme(
    
    legend.position = "none",               # Remove legend
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 20),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12, color = "black",face = "bold"),
    axis.title.x = element_text(size = 15, color = "black",face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing.x = unit(0.5, "cm")
  )

png("~/fig3_3_female.png",res=350,width=8000, height=4000)
fig3.3
dev.off()
######################################################################################################
df_plot <- bind_rows(df.all, df.men, df.women) %>%
  mutate(
    sex = factor(sex, levels = c("All", "Men", "Women")),
    # Ensure these exist / are consistent
    shape_group = as.factor(shape_group),
    is_pooled = as.factor(is_pooled),
    # nicer ordering within each score if you want
    study = fct_inorder(study),
    xmin_cap = pmax(lci, 0.3),
    xmax_cap = pmin(uci, 3),
    left_trunc  = lci < 0.3,
    right_trunc = uci > 3
  )

df_plot <- df_plot %>%
  filter(is.finite(or), is.finite(lci), is.finite(uci)) %>%
  filter(or > 0, lci > 0, uci > 0) %>%   # required for log10 scale
  droplevels()                            # remove empty score levels

df_plot_ss<-df_plot %>% filter(score %in% c("Acesulfame K", 
                                            "Aspartame", 
                                            "Saccharin", 
                                            "Sucralose"))

df_plot_sa<-df_plot %>% filter(score %in% c("Erythritol",
                                            "Maltitol",
                                            "Mannitol",
                                            "Sorbitol",
                                            "Xylitol"))

df_plot_ssg<-df_plot %>% filter(score %in% c("Sugar alcohols",
                                             "Synthetic sweeteners"))
################################################################################################################################
png("~/fig3_1_sens.png",res=350,width=8000, height=5000)
ggplot(
  df_plot_ss,
  aes(x = or, y = fct_rev(study), color = sex)
) +
  geom_point(aes(size = weight, shape = shape_group), position = position_dodge(width = 0.7), alpha = 0.9) +
  scale_size_continuous(range = c(1, 7))+
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(aes(xmin = lci, xmax = uci),
                 position = position_dodge(width = 0.7),
                 height = 0.15,
                 linewidth = 0.4) +
  guides(
    shape = "none",
    size  = "none",
    color = guide_legend(title = "Sex")
  )+
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  geom_text(
    data = df_plot_ss %>% filter(study == "Pooled"),
    aes(label = pooled_label, group = sex),
    position = position_dodge(width = 0.7), # 🔑 this aligns with points
    vjust = 0,         # vertical offset above diamond
    hjust= -0.3,
    size = 5,
    color = "black"
  )+
  scale_x_log10() +
  facet_wrap(~ score, scales = "free_y", ncol = 2, drop = TRUE) +
  theme_minimal() +
  theme(
    
    legend.position = "top",               # Remove legend
    legend.title = element_text(size = 24, face = "bold"),
    legend.text  = element_text(size = 22),
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 30),
    axis.text.y = element_text(size = 22, color = "black",face = "bold"),
    axis.text.x = element_text(size = 22, color = "black",face = "bold"),
    axis.title.x = element_text(size = 25, color = "black",face = "bold"),
    axis.title = element_text(size = 22),
    plot.title = element_text(face = "bold", hjust = 0.5) 
  )
dev.off()
################################################################################################################################
png("~/fig3_2_sens.png",res=350,width=10000, height=6000)
ggplot(
  df_plot_sa,
  aes(x = or, y = fct_rev(study), color = sex)
) +
  geom_point(aes(size = weight, shape = shape_group), position = position_dodge(width = 0.7), alpha = 0.9) +
  scale_size_continuous(range = c(1, 7))+
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(aes(xmin = lci, xmax = uci),
                 position = position_dodge(width = 0.7),
                 height = 0.15,
                 linewidth = 0.4) +
  guides(
    shape = "none",
    size  = "none",
    color = guide_legend(title = "Sex")
  )+
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  geom_text(
    data = df_plot_sa %>% filter(study == "Pooled"),
    aes(label = pooled_label, group = sex),
    position = position_dodge(width = 0.7), # 🔑 this aligns with points
    vjust = 0,         # vertical offset above diamond
    hjust= -0.2,
    size = 5,
    color = "black"
  )+
  scale_x_log10() +
  facet_wrap(~ score, scales = "free_y", ncol = 3, drop = TRUE) +
  theme_minimal() +
  theme(
    
    legend.position = "top",               # Remove legend
    legend.title = element_text(size = 24, face = "bold"),
    legend.text  = element_text(size = 22),
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 30),
    axis.text.y = element_text(size = 22, color = "black",face = "bold"),
    axis.text.x = element_text(size = 22, color = "black",face = "bold"),
    axis.title.x = element_text(size = 25, color = "black",face = "bold"),
    axis.title = element_text(size = 22),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
dev.off()
################################################################################################################################
png("~/fig3_3_sens.png",res=350,width=8000, height=4000)
ggplot(
  df_plot_ssg,
  aes(x = or, y = fct_rev(study), color = sex)
) +
  geom_point(aes(size = weight, shape = shape_group), position = position_dodge(width = 0.7), alpha = 0.9) +
  scale_size_continuous(range = c(1, 7))+
  scale_shape_manual(values = c("Square" = 15, "Diamond" = 18)) +  # square = study, diamond = pooled
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(aes(xmin = lci, xmax = uci),
                 position = position_dodge(width = 0.7),
                 height = 0.15,
                 linewidth = 0.4) +
  guides(
    shape = "none",
    size  = "none",
    color = guide_legend(title = "Sex")
  )+
  labs(
    x = "Odds Ratio (95% CI)",
    y = "",
    title = ""
  ) +
  geom_text(
    data = df_plot_ssg %>% filter(study == "Pooled"),
    aes(label = pooled_label, group = sex),
    position = position_dodge(width = 0.6), # 🔑 this aligns with points
    vjust = 0,         # vertical offset above diamond
    hjust= -0.2,
    size = 5,
    color = "black"
  )+
  scale_x_log10() +
  facet_wrap(~ score, scales = "free_y", ncol = 2, drop = TRUE) +
  theme_minimal() +
  theme(
    
    legend.position = "top",               # Remove legend
    legend.title = element_text(size = 24, face = "bold"),
    legend.text  = element_text(size = 22),
    panel.grid.major = element_blank(),     # Remove major gridlines
    panel.grid.minor = element_blank(),     # Remove minor gridlines
    panel.border = element_blank(),         # No border box
    axis.line = element_line(color = "black"),  # Keep x and y axis lines
    axis.ticks = element_line(color = "black"), # Keep tick marks
    axis.text = element_text(color = "black",face = "bold"),  # Keep axis labels
    strip.text = element_text(face = "bold",size = 30),
    axis.text.y = element_text(size = 22, color = "black",face = "bold"),
    axis.text.x = element_text(size = 22, color = "black",face = "bold"),
    axis.title.x = element_text(size = 25, color = "black",face = "bold"),
    axis.title = element_text(size = 22),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
dev.off()
################################################################################################################################
sample_size<-read.csv("~/replicate_bb4_score_sex_cc.csv",header=T)

table_cc <- sample_size %>%
  mutate(
    sex = recode(sex, Male = "Men", Female = "Women"),
    case_control = recode(case_control, Con = "Control", Case = "T2D")
  ) %>%
  # reshape to wide: Control and T2D as columns
  pivot_wider(
    names_from = case_control,
    values_from = freq,
    values_fill = 0
  ) %>%
  # create formatted strings
  mutate(
    case_control = paste0(T2D, "/", Control)
  ) %>%
  select(study, sex, case_control) %>%
  # reshape again: Men/Women as columns
  pivot_wider(
    names_from = sex,
    values_from = case_control,
    names_prefix = "case_control_"
  ) %>%
  arrange(study)

write.csv(table_cc,"~/table_cc.csv",row.names=F)
################################################################################################################################
## Address heterogeneity ##
################################################################################################################################
df.all.hetero.1 <- df.all %>%
  filter(!study %in% "Pooled") %>%
  group_by(score) %>%
  summarise(
    # Perform meta-analysis using the current group's data
    meta = list(rma(yi = beta, sei = se, data = cur_data(), method = "REML")),
    # Extract and calculate metrics
    pooled_beta = meta[[1]]$b[1],
    pooled_se = meta[[1]]$se,
    pooled_or = exp(pooled_beta),
    pooled_lci = exp(pooled_beta - 1.96 * pooled_se),
    pooled_uci = exp(pooled_beta + 1.96 * pooled_se),
    pooled_p = meta[[1]]$pval,
    tau2 = meta[[1]]$tau2,
    i2 = meta[[1]]$I2,
    q = meta[[1]]$QE,
    q_p = meta[[1]]$QEp,
    .groups = 'drop'
  ) %>%
  select(score,i2,q,q_p)  # Remove the list column

df.all.hetero.2 <- df.all %>%
  filter(!study %in% "Pooled") %>%
  group_by(score) %>%
  summarise(
    # Perform meta-analysis using the current group's data
    meta = list(rma(yi = beta, sei = se, data = cur_data(), method = "FE")),
    # Extract and calculate metrics
    pooled_beta = meta[[1]]$b[1],
    pooled_se = meta[[1]]$se,
    pooled_or = exp(pooled_beta),
    pooled_lci = exp(pooled_beta - 1.96 * pooled_se),
    pooled_uci = exp(pooled_beta + 1.96 * pooled_se),
    pooled_p = meta[[1]]$pval,
    tau2 = meta[[1]]$tau2,
    i2 = meta[[1]]$I2,
    q = meta[[1]]$QE,
    q_p = meta[[1]]$QEp,
    .groups = 'drop'
  ) %>%
  mutate(OR_CI=sprintf("%.2f(%.2f,%.2f)", pooled_or, pooled_lci, pooled_uci),
         sex="All") %>%
  select(score,sex,OR_CI)  # Remove the list column

df.all.hetero <- df.all.hetero.2 %>%
  left_join(df.all.hetero.1,by="score")


df.men.hetero.1 <- df.men %>%
  filter(!study %in% "Pooled") %>%
  group_by(score) %>%
  summarise(
    # Perform meta-analysis using the current group's data
    meta = list(rma(yi = beta, sei = se, data = cur_data(), method = "REML")),
    # Extract and calculate metrics
    pooled_beta = meta[[1]]$b[1],
    pooled_se = meta[[1]]$se,
    pooled_or = exp(pooled_beta),
    pooled_lci = exp(pooled_beta - 1.96 * pooled_se),
    pooled_uci = exp(pooled_beta + 1.96 * pooled_se),
    pooled_p = meta[[1]]$pval,
    tau2 = meta[[1]]$tau2,
    i2 = meta[[1]]$I2,
    q = meta[[1]]$QE,
    q_p = meta[[1]]$QEp,
    .groups = 'drop'
  ) %>%
  select(score,i2,q,q_p)  # Remove the list column

df.men.hetero.2 <- df.men %>%
  filter(!study %in% "Pooled") %>%
  group_by(score) %>%
  summarise(
    # Perform meta-analysis using the current group's data
    meta = list(rma(yi = beta, sei = se, data = cur_data(), method = "FE")),
    # Extract and calculate metrics
    pooled_beta = meta[[1]]$b[1],
    pooled_se = meta[[1]]$se,
    pooled_or = exp(pooled_beta),
    pooled_lci = exp(pooled_beta - 1.96 * pooled_se),
    pooled_uci = exp(pooled_beta + 1.96 * pooled_se),
    pooled_p = meta[[1]]$pval,
    tau2 = meta[[1]]$tau2,
    i2 = meta[[1]]$I2,
    q = meta[[1]]$QE,
    q_p = meta[[1]]$QEp,
    .groups = 'drop'
  ) %>%
  mutate(OR_CI=sprintf("%.2f(%.2f,%.2f)", pooled_or, pooled_lci, pooled_uci),
         sex="Men") %>%
  select(score,sex,OR_CI)  # Remove the list column

df.men.hetero <- df.men.hetero.2 %>%
  left_join(df.men.hetero.1,by="score")

df.women.hetero.1 <- df.women %>%
  filter(!study %in% "Pooled") %>%
  group_by(score) %>%
  summarise(
    # Perform meta-analysis using the current group's data
    meta = list(rma(yi = beta, sei = se, data = cur_data(), method = "REML")),
    # Extract and calculate metrics
    pooled_beta = meta[[1]]$b[1],
    pooled_se = meta[[1]]$se,
    pooled_or = exp(pooled_beta),
    pooled_lci = exp(pooled_beta - 1.96 * pooled_se),
    pooled_uci = exp(pooled_beta + 1.96 * pooled_se),
    pooled_p = meta[[1]]$pval,
    tau2 = meta[[1]]$tau2,
    i2 = meta[[1]]$I2,
    q = meta[[1]]$QE,
    q_p = meta[[1]]$QEp,
    .groups = 'drop'
  ) %>%
  select(score,i2,q,q_p)  # Remove the list column

df.women.hetero.2 <- df.women %>%
  filter(!study %in% "Pooled") %>%
  group_by(score) %>%
  summarise(
    # Perform meta-analysis using the current group's data
    meta = list(rma(yi = beta, sei = se, data = cur_data(), method = "FE")),
    # Extract and calculate metrics
    pooled_beta = meta[[1]]$b[1],
    pooled_se = meta[[1]]$se,
    pooled_or = exp(pooled_beta),
    pooled_lci = exp(pooled_beta - 1.96 * pooled_se),
    pooled_uci = exp(pooled_beta + 1.96 * pooled_se),
    pooled_p = meta[[1]]$pval,
    tau2 = meta[[1]]$tau2,
    i2 = meta[[1]]$I2,
    q = meta[[1]]$QE,
    q_p = meta[[1]]$QEp,
    .groups = 'drop'
  ) %>%
  mutate(OR_CI=sprintf("%.2f(%.2f,%.2f)", pooled_or, pooled_lci, pooled_uci),
         sex="Women") %>%
  select(score,sex, OR_CI)  # Remove the list column

df.women.hetero <- df.women.hetero.2 %>%
  left_join(df.women.hetero.1,by="score")


df.hetero <- 
  df.all.hetero %>%
  left_join(df.men.hetero,by="score") %>%
  left_join(df.women.hetero,by="score") 

write.csv(df.hetero,"~/hetero.csv",row.names=F)
################################################################################################################################
## Address confounding ##
################################################################################################################################
df.unadjust <- read.csv("~/replicate_bb4_score_unadjut.csv",header=T)

df.unadjust$study<-ifelse(df.unadjust$study=="pooled","Pooled",df.unadjust$study)

# Create flag for truncated CIs
df.unadjust <- df.unadjust %>%
  mutate(
    OR_CL_Unadjust = ifelse(
      study == "Pooled",
      sprintf("%.2f(%.2f,%.2f)", or, lci, uci),
      NA
    )
  )%>%
  filter(study=="Pooled")%>%
  select(score,OR_CL_Unadjust)

df.nobmi <- read.csv("~/replicate_bb4_score_nobmi.csv",header=T)

df.nobmi$study<-ifelse(df.nobmi$study=="pooled","Pooled",df.nobmi$study)

# Create flag for truncated CIs
df.nobmi <- df.nobmi %>%
  mutate(
    OR_CL_nobmi = ifelse(
      study == "Pooled",
      sprintf("%.2f(%.2f,%.2f)", or, lci, uci),
      NA
    )
  )%>%
  filter(study=="Pooled")%>%
  select(score,OR_CL_nobmi)

df.full<-df.all.hetero %>%
  select(score,OR_CI)

df.3sens<-
  df.unadjust %>%
  left_join(df.nobmi,by="score") %>%
  left_join(df.full,by="score")

write.csv(df.3sens,"~/replciate_3sens.csv",row.names = F)


