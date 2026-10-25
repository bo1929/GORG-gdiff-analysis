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
library(readr)
library(cowplot); library(patchwork); library(stringr)
library(ggstar);
library(grid)

dfg <- vroom("results/genes-joined-all.tsv.gz")
dff <- dfg
dfs <- vroom("results/genome-all.tsv.gz")
dfr <- vroom("resource-scaling-nt8-10x100.tsv")

p1 <- dfr %>%
  filter(case %in% c("NxQ", "index_N")) %>%
  filter(threads == 8) %>%
  mutate(case = case_when(
    case == "NxQ" ~ "Query\n10x100",
    case == "index_N" ~ "Index\n100",
    TRUE ~ case
  )) %>%
  ggplot(aes(y = case, x = wall_s, fill = reorder(tool, -wall_s), group = reorder(tool, -wall_s))) +
  geom_col(width=0.5, position = position_dodge(width = 0.66), aes(color=reorder(tool, -wall_s))) +
  geom_text(
    aes(label = round(wall_s, 0), color=reorder(tool, -wall_s)),
    position = position_dodge2(width = 0.66),
    hjust = -0.15,
    size = 3,
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
    case == "NxQ" ~ "Query\n10x100",
    case == "index_N" ~ "Index\n100",
    TRUE ~ case
  )) %>%
  ggplot(aes(y = case, x = peak_rss_mb/1000, fill = reorder(tool, -wall_s), group = reorder(tool, -wall_s))) +
  geom_col(width=0.5, position = position_dodge(width = 0.66), aes(color=reorder(tool, -wall_s))) +
  geom_text(
    aes(label = round(peak_rss_mb/1000, 3), color=reorder(tool, -wall_s)),
    position = position_dodge(width = 0.66),
    size=3,
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

dfg <- dff
thw <- 0.5

inter <- dfg %>%
  filter(included == 1, wf_cov >= thw, is.finite(wf_rate), is.finite(d_gd))

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


options(scipen = 999)

df <- vroom("threshold-per-pair.tsv", comment = "#")
df <- df %>% mutate(ratevar=if_else(grepl("a5", variant), "High rate variation", "Low rate variation"))
df <- df %>% filter(dataset == "uniform-intergenic")

pth <- df %>% 
  filter(config=="k25-w29-frac1.0") %>%
  filter(!grepl("gnd00000", variant)) %>%
  filter(l == 1000) %>%
  ggplot() +
  # geom_point() +
  facet_wrap(~ratevar, nrow = 2) +
  stat_summary_bin(size=0.15, alpha=0.75, aes(y = threshold_low, x = true_gnd, color=factor(level)), bins = 20) +
  stat_summary_bin(size=0.15, alpha=0.75, aes(y = threshold_high, x = true_gnd, color=factor(level)), bins = 20) +
  geom_abline(color="red") +
  scale_color_discrete(palette="viridis") +
  labs(x=TeX(r'($D$)'), y=TeX(r'(Threshold at $\alpha$)'), color=TeX(r'($\alpha$)')) +
  theme_bw()
pth

plot_grid(
  p3 + theme(legend.position = "bottom"),
  p1 + theme(legend.position = "none", axis.title.y = element_blank()),
  p2 + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank()),
  pth + theme(legend.position = "bottom") + 
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.margin = margin(0, 0, 0, 0)),
  nrow=1, rel_widths = c(0.8, 0.75, 0.55, 0.7)
)
ggsave2("./combined-plot.pdf", width = 9, height = 3.5)

draw_key_arrow <- function(data, params, size) {
  grid::gTree(children = grid::gList(
    grid::textGrob("empirical", x = 0.00, hjust = 0.25, gp = grid::gpar(fontsize = 9)),
    grid::segmentsGrob(0.38, 0.5, 0.62, 0.5,
                       gp = grid::gpar(lwd = 1.2),
                       arrow = grid::arrow(length = grid::unit(2, "pt"), type = "closed")),
    grid::textGrob("model", x = 0.98, hjust = 1, gp = grid::gpar(fontsize = 9))
  ))
}

arrow_legend <- data.frame(
  ratevar = unique(df$ratevar)[1],
  x = 0, xend = 0, y = 0, yend = 0,
  direction = "arrow"
)

p2 <- df %>%
  filter(config=="k25-w29-frac1.0") %>%
  filter(!grepl("gnd00000", variant)) %>%
  filter(l == 1000) %>%
  # mutate(gnd_bin=cut(true_gnd, c(1:20)/100)) %>%
  mutate(gnd_bin=round(true_gnd, 2)) %>%
  filter(level==0.01) %>%
  group_by(ratevar, gnd_bin, level) %>%
  summarize(
    threshold_low=mean(threshold_low),
    threshold_high=mean(threshold_high),
    truth_realized_tail_low=mean(truth_realized_tail_low),
    truth_realized_tail_high=mean(truth_realized_tail_high)
  ) %>%
  ggplot() +
  # geom_abline() +
  # geom_point(alpha=0.1) +
  # facet_grid(rows=vars(ratevar), cols=vars(level)) +\
  facet_wrap(~ratevar, scale="free") +
  geom_segment(arrow = arrow(length = unit(2, "pt"), type = "closed"),
               aes(color=gnd_bin, xend = threshold_low, x = truth_realized_tail_low, y=threshold_high, yend = truth_realized_tail_high)
               # aes(color=gnd_bin, x = threshold_low, y = truth_realized_tail_low, xend=threshold_high, yend = truth_realized_tail_high)
  ) +
  # legend-only arrow; colour=NA keeps it out of the panel
  geom_segment(data=arrow_legend,
               aes(x=x, xend=xend, y=y, yend=yend, linetype=direction),
               inherit.aes=FALSE, colour=NA,
               arrow=arrow(length=unit(2, "pt"), type="closed"),
               key_glyph=draw_key_arrow, show.legend=TRUE) +
  scale_linetype_manual(values=c("arrow"="solid"), name=NULL) +
  guides(linetype=guide_legend(label=FALSE, title=NULL, keywidth=unit(90, "pt"), keyheight=unit(0.5, "cm"), order=2)) +
  scale_color_viridis_c(option = "mako", end = 0.85) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x=TeX(r'(Lower tail ($\alpha=0.01$))'), y=TeX(r'(Upper tail ($\alpha=0.01$))'), color=TeX(r'($D$)')) # +
# theme(strip.text = element_blank(), strip.background = element_blank())
p2
