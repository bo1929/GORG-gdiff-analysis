library(dplyr); library(ggplot2); library(vroom); library(latex2exp)
library(scales)
library(ggpubr)
library(cowplot)
library(tidyr)
library(ggbreak) 

# library(ggforce)
# install.packages("ggforce")
# https://ggforce.data-imaginist.com

mc = c( "#9D4030", "#D0D55C", "#C3A97E", "#809D6F", "#769DA6", "#A69E33", "#734002", "#595622", "#8C873F",  "#474B71")
lvl=c("gdiff", "fastani", "mash", "dashing2", "skani")

df11 <- vroom("output-11/resources.tsv") %>% mutate(seed = 11)
df17 <- vroom("output-17/resources.tsv") %>% mutate(seed = 17)
df19 <- vroom("output-19/resources.tsv") %>% mutate(seed = 19)
df29 <- vroom("output-29/resources.tsv") %>% mutate(seed = 29)
df79 <- vroom("output-79/resources.tsv") %>% mutate(seed = 79)
df <- rbind(df11, df17, df19, df29, df79)

df <- df %>% mutate(wall_sec=if_else(method=="skani" & phase=="sketch", 0, wall_sec)) %>%
  group_by(method, seed) %>%
  summarise(wall_sec=sum(wall_sec, na.rm=T), maxrss_mib=max(maxrss_mib, na.rm=T))

df <- df %>% mutate(method=factor(method, levels=lvl))
df %>% 
  # complete(phase, method, fill = list(val = NA)) %>%
  mutate(lbl=if_else(is.na(wall_sec), "x", "")) %>%
  mutate(wall_sec=if_else(is.na(wall_sec), 0, wall_sec)) %>%
  ggplot() +
  aes(x=method, y=wall_sec, fill=method) +
  # facet_wrap(~phase) +
  stat_summary(geom="bar", show.legend = F) +
  geom_text(aes(label=lbl, y=0.5), size=4, color="black") +
  theme_minimal_grid(font_size = 10) +
  # scale_y_break(c(10, 137.5), symbol = "slash") +
  scale_fill_manual(values=mc) +
  # coord_cartesian(ylim=c(0, 142)) +
  labs(x="Method", y="Wall time (sec)") +
  # scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 138, 140)) +
  # scale_y_continuous(transform = "log1p", breaks = c(1, 5, 15, 30, 60, 120)) +
  scale_y_continuous(transform = "sqrt", breaks = c(1, 10, 30, 60, 120)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))
ggsave2("./R-running_time-ani.pdf", width = 3, height = 3)

rt <- df %>% 
  # complete(phase, method, fill = list(val = NA)) %>%
  mutate(lbl=if_else(is.na(wall_sec), "x", "")) %>%
  mutate(wall_sec=if_else(is.na(wall_sec), 0, wall_sec)) %>%
  ggplot() +
  aes(y=method, x=wall_sec, fill=method) +
  # facet_wrap(~phase) +
  stat_summary(geom="bar", show.legend = F) +
  geom_text(aes(label=lbl, y=0.5), size=4, color="black") +
  theme_minimal() +
  # scale_y_break(c(10, 137.5), symbol = "slash") +
  scale_fill_manual(values=mc) +
  # coord_cartesian(ylim=c(0, 142)) +
  labs(y="", x="Wall time (sec)") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  # scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 138, 140)) +
  scale_x_continuous(transform = "sqrt", breaks = c(1, 5, 15, 45, 90)) # +
  # theme(axis.text.x = element_text(angle=45, hjust = 1))
rt
ggsave2("./R-running_time-ani.pdf", width = 2.25, height = 2)


df %>% 
  # complete(phase, method, fill = list(val = NA)) %>%
  mutate(lbl=if_else(is.na(wall_sec), "x", "")) %>%
  mutate(wall_sec=if_else(is.na(wall_sec), 0, wall_sec)) %>%
  ggplot() +
  aes(x=reorder(method, -wall_sec), y=maxrss_mib, fill=method) +
  # facet_wrap(~phase) +
  stat_summary(geom="bar", show.legend = F) +
  geom_text(aes(label=lbl, y=75), size=4, color="black") +
  theme_minimal_grid(font_size = 10) +
  # scale_y_break(c(7.5, 137.5), symbol = "slash") +
  scale_fill_manual(values=mc) +
  # coord_cartesian(ylim=c(0, 142)) +
  # scale_y_continuous(transform = "log1p", breaks = c(5, 50, 250, 1000)) +
  labs(x="Method", y="Peak memory (MB)") # +
  # theme(axis.text.x = element_text(angle=45, hjust = 1))
ggsave2("./R-peak_memory-ani.pdf", width = 3, height = 3)

mem <- df %>% 
  # complete(phase, method, fill = list(val = NA)) %>%
  mutate(lbl=if_else(is.na(wall_sec), "x", "")) %>%
  mutate(wall_sec=if_else(is.na(wall_sec), 0, wall_sec)) %>%
  ggplot() +
  aes(y=method, x=maxrss_mib, fill=method) +
  # facet_wrap(~phase) +
  stat_summary(geom="bar", show.legend = F) +
  # geom_text(aes(label=lbl, x=75), size=4, color="black") +
  theme_minimal() +
  # scale_y_break(c(7.5, 137.5), symbol = "slash") +
  scale_fill_manual(values=mc) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  # coord_cartesian(ylim=c(0, 142)) +
  # scale_y_continuous(transform = "log1p", breaks = c(5, 50, 250, 1000)) +
  labs(y="", x="Peak memory (MB)") #+
  # theme(axis.text.x = element_text(angle=45, hjust = 1))
mem
ggsave2("./R-peak_memory-ani.pdf", width = 2.25, height = 2)

plot_grid(rt, mem + theme(axis.text.y = element_blank(), axis.title.y = element_blank()))
ggsave2("./resource-use-ani.pdf", width = 4, height = 1.75)
