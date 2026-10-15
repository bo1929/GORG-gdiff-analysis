library(vroom)
library(rstatix)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(latex2exp)
library(cowplot)
# library(fitdistrplus)
# install.packages("ggforce")
library(ggforce)

dfg <- vroom("results/genes-joined-all.tsv.gz")
dfs <- vroom("results/genome-all.tsv.gz")
dfr <- vroom("resource-scaling-nt8-10x100.tsv")

dfr %>%
  filter(case %in% c("NxQ", "index_N")) %>%
  filter(threads == 8) %>%
  mutate(case = case_when(
    case == "NxQ" ~ "Map 10x100",
    case == "index_N" ~ "Index 100",
    TRUE ~ case
  )) %>%
  ggplot(aes(x = case, y = wall_s, fill = reorder(tool, -wall_s), group = reorder(tool, -wall_s))) +
  geom_col(width=0.5, position = position_dodge(width = 0.66), aes(color=reorder(tool, -wall_s))) +
  geom_text(
    aes(label = round(wall_s, 1), color=reorder(tool, -wall_s)), 
    position = position_dodge(width = 0.66), 
    vjust = -0.5, 
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#595622", "#9D4030")) +
  scale_color_manual(values = c("#595622", "#9D4030")) +
  theme_minimal_hgrid() +
  scale_y_continuous(n.breaks = 6) +
  labs(x = "Stage", y = "Running time (sec)", fill = "Method", color="Method") # +
  # coord_cartesian(ylim = c(0, 60))
ggsave2("./S-running-time-wfmash-nt8.pdf", width = 4, height = 4.5)

dfr %>%
  filter(case %in% c("NxQ", "index_N")) %>%
  filter(threads == 8) %>%
  mutate(case = case_when(
    case == "NxQ" ~ "Map 10x100",
    case == "index_N" ~ "Index 100",
    TRUE ~ case
  )) %>%
  ggplot(aes(x = case, y = peak_rss_mb, fill = reorder(tool, -wall_s), group = reorder(tool, -wall_s))) +
  geom_col(width=0.5, position = position_dodge(width = 0.66), aes(color=reorder(tool, -wall_s))) +
  geom_text(
    aes(label = round(peak_rss_mb, 1), color=reorder(tool, -wall_s)), 
    position = position_dodge(width = 0.66), 
    size=4,
    vjust = -0.5, 
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#595622", "#9D4030")) +
  scale_color_manual(values = c("#595622", "#9D4030")) +
  theme_minimal_hgrid() +
  scale_y_continuous(n.breaks = 6) +
  labs(x = "Stage", y = "Peak memory (MB)", fill = "Method", color="Method")
ggsave2("./S-peak-memory-wfmash-nt8.pdf", width = 4.25, height = 4.5)

thw <- 0.0
ggplot() +
  geom_abline() +
  stat_cor(
    data = dfg %>%
      filter(included == 1, wf_cov > thw), aes(label = ..r.label.., x=truth, y=wf_rate), color="#595622", label.y = 0.36) +
  stat_summary_bin(
    data = dfg %>%
      filter(included == 1, wf_cov > thw), aes(x=truth, y=wf_rate), alpha=0.85, color="#595622") +
  stat_cor(
    data = dfg %>%
      filter(included == 1), aes(label = ..r.label.., x=truth, y=d_gd), color="#9D4030", label.y = 0.33) +
  stat_summary_bin(
    data = dfg %>%
      filter(included == 1), aes(x=truth, y=d_gd), alpha=0.85, shape=17, color="#9D4030") +
  geom_text(
    data = function(d) dfg %>% filter(included == 1, wf_cov > thw) %>%
      summarize(mape = mean(abs((wf_rate - truth) / (truth))) * 100,  .groups = "drop") %>% 
      mutate(x = 0.075, y = 0.355),
    aes(x=x, y=y, label=paste0("MAPE = ", round(mape, 2), "%")),
    hjust=0, vjust=0, show.legend = F, color="#595622"
  ) + 
  geom_text(
    data = function(d) dfg %>% filter(included == 1) %>%
      summarize(mape = mean(abs((d_gd - truth) / (truth))) * 100,  .groups = "drop") %>% 
      mutate(x = 0.075, y = 0.325),
    aes(x=x, y=y, label=paste0("MAPE = ", round(mape, 2), "%")),
    hjust=0, vjust=0, show.legend = F, color="#9D4030"
  ) +
  theme_bw() +
  coord_cartesian(xlim=c(0, 0.4), ylim=c(0, 0.4)) +
  labs(x=TeX(r'(Gene $GND$)'), y=TeX(r'($\hat{GND}$)'))
ggsave2("./S-mape-corr-wfmash.pdf", width = 4.25, height = 4.5)

dfg %>% 
  filter(included == 1) %>% 
  filter(truth < 0.5) %>%
  mutate(gene_gnd = cut(truth, c(0, 0.05, 0.1, 0.2, 0.3, 0.5))) %>% 
  group_by(gene_gnd) %>% 
  # Calculate proportions and MAPE for each bin
  summarise(
    nwfmash = sum(included == 1 & wf_cov > 0) / n(), 
    ngdiff  = sum(included == 1 & is.finite(d_gd)) / n(),
    mape_wf = mean(abs((wf_rate[wf_cov > 0.0] - truth[wf_cov > 0.0]) / truth[wf_cov > 0.0]), na.rm = TRUE),
    mape_gd = mean(abs((d_gd - truth) / truth), na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  ggplot(aes(x = gene_gnd)) +
  geom_col(aes(y = mape_wf), fill = "#595622", color="black", alpha = 0.99, width = 0.33, position = position_nudge(x = -0.22)) +
  geom_col(aes(y = mape_gd), fill = "#9D4030", color="black", alpha = 0.99, width = 0.33, position = position_nudge(x = 0.22)) +
  geom_text(
    aes(y = mape_wf, label = paste0(round(nwfmash*100), "%")), 
    vjust = -0.5, position = position_nudge(x = -0.22), color = "#595622", size = 3, fontface = "bold"
  ) +
  geom_text(
    aes(y = mape_gd, label = paste0(round(ngdiff*100), "%")), 
    vjust = -0.5, position = position_nudge(x = 0.22), color = "#9D4030", size = 3, fontface = "bold"
  ) +
  theme_minimal_hgrid() +
  scale_y_continuous(labels = percent) +
  labs(x = "Gene distances", y = "MAPE")
ggsave2("./S-mape-cov-wfmash.pdf", width = 4.5, height = 4.5)
