library(vroom)
library(rstatix)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(scales)
library(latex2exp)
library(cowplot)
# library(fitdistrplus)
# install.packages("ggforce")
library(ggforce)
library(ggpubr)

dfg <- vroom("results/genes-joined-all.tsv.gz")
dff <- dfg
dfs <- vroom("results/genome-all.tsv.gz")
dfr <- vroom("resource-scaling-nt8-10x100.tsv")

p1 <- dfr %>%
  filter(case %in% c("NxQ", "index_N")) %>%
  filter(threads == 8) %>%
  mutate(case = case_when(
    case == "NxQ" ~ "Query 10x100",
    case == "index_N" ~ "Index 100",
    TRUE ~ case
  )) %>%
  ggplot(aes(y = case, x = wall_s, fill = reorder(tool, -wall_s), group = reorder(tool, -wall_s))) +
  geom_col(width=0.5, position = position_dodge(width = 0.66), aes(color=reorder(tool, -wall_s))) +
  geom_text(
    aes(label = round(wall_s, 0), color=reorder(tool, -wall_s)),
    position = position_dodge2(width = 0.66),
    hjust = -0.15,
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#595622", "#9D4030")) +
  scale_color_manual(values = c("#595622", "#9D4030")) +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color="lightgray")) +
  scale_x_continuous(n.breaks = 6) +
  labs(y = "Stage", x = "Running time (sec)", fill = "Method", color="Method") +
  coord_cartesian(xlim = c(0, 1800))
p1
ggsave2("./S-running-time-wfmash-nt8.pdf", width = 4, height = 4.5)

p2 <- dfr %>%
  filter(case %in% c("NxQ", "index_N")) %>%
  filter(threads == 8) %>%
  mutate(case = case_when(
    case == "NxQ" ~ "Query 10x100",
    case == "index_N" ~ "Index 100",
    TRUE ~ case
  )) %>%
  ggplot(aes(y = case, x = peak_rss_mb/1000, fill = reorder(tool, -wall_s), group = reorder(tool, -wall_s))) +
  geom_col(width=0.5, position = position_dodge(width = 0.66), aes(color=reorder(tool, -wall_s))) +
  geom_text(
    aes(label = round(peak_rss_mb/1000, 3), color=reorder(tool, -wall_s)),
    position = position_dodge(width = 0.66),
    size=4,
    hjust = -0.15,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#595622", "#9D4030")) +
  scale_color_manual(values = c("#595622", "#9D4030")) +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color="lightgray")) +
  scale_x_continuous(n.breaks = 6) +
  labs(y = "Stage", x = "Peak memory (GB)", fill = "Method", color="Method") +
  coord_cartesian(xlim = c(0, 5.5))
p2
ggsave2("./S-peak-memory-wfmash-nt8.pdf", width = 4.25, height = 4.5)

# no truth filter: MAE needs no division, so genes at D = 0 stay in.
dfg <- dff
thw <- 0.5

# Both methods are scored on the same genes: wfmash aligned >= 50% of the gene AND
# gdiff has an estimate (the intersection).
inter <- dfg %>%
  filter(included == 1, wf_cov >= thw, is.finite(wf_rate), is.finite(d_gd))

# Coverage denominator: everything except the length filter and the no-window-inside case.
# The rest counts as a miss for a method that has no estimate: gdiff is uncovered when d_gd
# is NA (gdiff_nan), wfmash is uncovered when wf_cov < 0.5 (including no alignment at all).
eligible <- dfg %>% filter(!reason %in% c("shorter_than_window", "no_window_inside"))
cov_wf <- mean(eligible$wf_cov >= thw & is.finite(eligible$wf_rate))
cov_gd <- mean(is.finite(eligible$d_gd))

mae_wf <- mean(abs(inter$wf_rate - inter$truth))
mae_gd <- mean(abs(inter$d_gd - inter$truth))
r_wf <- cor(inter$truth, inter$wf_rate)
r_gd <- cor(inter$truth, inter$d_gd)
lim <- 1.05 * max(c(inter$truth, inter$wf_rate, inter$d_gd), na.rm = TRUE)

p3 <- ggplot() +
  geom_abline() +
  stat_summary_bin(
    data = inter, aes(x=truth, y=wf_rate, colour="wfmash"), alpha=0.85) +
  stat_summary_bin(
    data = inter, aes(x=truth, y=d_gd, colour="gdiff"), alpha=0.85, shape=17) +
  annotate("text", x=0.02*lim, y=0.90*lim, hjust=0, colour="#595622", size=3.1,
           label=sprintf("wfmash  R = %.4f\ncoverage = %.1f%%", r_wf, 100*cov_wf)) +
  annotate("text", x=0.02*lim, y=0.75*lim, hjust=0, colour="#9D4030", size=3.1,
           label=sprintf("gdiff   R = %.4f\ncoverage = %.1f%%", r_gd, 100*cov_gd)) +
  scale_colour_manual(values = c(wfmash="#595622", gdiff="#9D4030"), name = NULL) +
  theme_bw() +
  coord_cartesian(xlim=c(0, lim), ylim=c(0, lim)) +
  labs(x=TeX(r'(Gene $D$)'), y=TeX(r'(Gene $\hat{D}$)')) +
  coord_cartesian(xlim=c(0, 0.41), ylim=c(0, 0.41))
p3
ggsave2("./S-mape-corr-wfmash.pdf", width = 4.25, height = 4.5)

# MAE per gene-D bin on the same intersection; the labels above the bars are the
# share of eligible genes in that bin for which the method has an estimate.
brks <- c(0, 0.05, 0.1, 0.2, 0.3, 0.55)

cov_bin <- eligible %>%
  mutate(gene_gnd = cut(truth, brks, include.lowest = TRUE)) %>%
  group_by(gene_gnd) %>%
  summarise(
    n_all = n(),
    nwf = sum(wf_cov >= thw & is.finite(wf_rate)),
    ngd = sum(is.finite(d_gd)),
    .groups = "drop"
  )

mae_bin <- inter %>%
  mutate(gene_gnd = cut(truth, brks, include.lowest = TRUE)) %>%
  group_by(gene_gnd) %>%
  summarise(
    mae_wf = mean(abs(wf_rate - truth)),
    mae_gd = mean(abs(d_gd - truth)),
    n = n(),
    .groups = "drop"
  )

bins <- left_join(mae_bin, cov_bin, by = "gene_gnd")

bins %>%
  ggplot(aes(x = gene_gnd)) +
  geom_col(aes(y = mae_wf), fill = "#595622", color="black", alpha = 0.99, width = 0.33, position = position_nudge(x = -0.22)) +
  geom_col(aes(y = mae_gd), fill = "#9D4030", color="black", alpha = 0.99, width = 0.33, position = position_nudge(x = 0.22)) +
  geom_text(
    aes(y = mae_wf, label = paste0(round(100 * nwf / n_all), "%")),
    vjust = -0.5, position = position_nudge(x = -0.22), color = "#595622", size = 3, fontface = "bold"
  ) +
  geom_text(
    aes(y = mae_gd, label = paste0(round(100 * ngd / n_all), "%")),
    vjust = -0.5, position = position_nudge(x = 0.22), color = "#9D4030", size = 3, fontface = "bold"
  ) +
  theme_minimal_hgrid() +
  scale_y_continuous(labels = percent) +
  labs(x = TeX(r'(Gene $D$)'), y = "MAE")
ggsave2("./S-mape-cov-wfmash.pdf", width = 4.5, height = 4.5)

plot_grid(
  p3 + theme(legend.position = "bottom"),
  p1 + theme(legend.position = "none", axis.title.y = element_blank()),
  p2 + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank()),
  nrow=1, rel_widths = c(0.8, 0.75, 0.5)
)
