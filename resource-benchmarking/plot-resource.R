library(dplyr); library(ggplot2); library(vroom); library(latex2exp)
library(scales)
library(ggpubr)
library(cowplot)
library(tidyr)
library(ggbreak) 

# library(ggforce)
# install.packages("ggforce")
# https://ggforce.data-imaginist.com

mc = c("#809D6F", "#D0D55C", "#9D4030", "#C3A97E", "#769DA6", "#A69E33", "#734002", "#595622", "#8C873F",  "#474B71")

df11 <- vroom("output-11/resources.tsv") %>% mutate(seed = 11)
df17 <- vroom("output-17/resources.tsv") %>% mutate(seed = 17)
df19 <- vroom("output-19/resources.tsv") %>% mutate(seed = 19)
df29 <- vroom("output-29/resources.tsv") %>% mutate(seed = 29)
df79 <- vroom("output-79/resources.tsv") %>% mutate(seed = 79)
df <- rbind(df11, df17, df19, df29, df79)

df %>% 
  complete(phase, method, fill = list(val = NA)) %>%
  mutate(lbl=if_else(is.na(wall_sec), "x", "")) %>%
  mutate(wall_sec=if_else(is.na(wall_sec), 0, wall_sec)) %>%
  ggplot() +
  aes(x=method, y=wall_sec, fill=method) +
  facet_wrap(~phase) +
  stat_summary(geom="bar", show.legend = F) +
  geom_text(aes(label=lbl, y=0.5), size=4, color="black") +
  theme_minimal_grid(font_size = 10) +
  scale_y_break(c(7.5, 137.5), symbol = "slash") +
  scale_fill_manual(values=mc) +
  coord_cartesian(ylim=c(0, 142)) +
  labs(x="Method", y="Wall time (sec)") +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 138, 140)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))
ggsave2("./R-running_time-ani.pdf", width = 3, height = 3)

df %>% 
  complete(phase, method, fill = list(val = NA)) %>%
  mutate(lbl=if_else(is.na(wall_sec), "x", "")) %>%
  mutate(wall_sec=if_else(is.na(wall_sec), 0, wall_sec)) %>%
  ggplot() +
  aes(x=method, y=maxrss_mib, fill=method) +
  facet_wrap(~phase) +
  stat_summary(geom="bar", show.legend = F) +
  geom_text(aes(label=lbl, y=75), size=4, color="black") +
  theme_minimal_grid(font_size = 10) +
  # scale_y_break(c(7.5, 137.5), symbol = "slash") +
  scale_fill_manual(values=mc) +
  # coord_cartesian(ylim=c(0, 142)) +
  labs(x="Method", y="Peak memory (MB)") +
  theme(axis.text.x = element_text(angle=45, hjust = 1))
ggsave2("./R-peak_memory-ani.pdf", width = 3, height = 3)
