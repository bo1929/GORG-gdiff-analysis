library(dplyr); library(ggplot2); library(vroom); library(latex2exp)
library(scales)
library(ggpubr)
library(cowplot)
lvl=c("gdiff", "fastani", "mash", "dashing2", "skani")
mc = c("#809D6F", "#D0D55C", "#BD4030", "#C3A97E", "#769DA6", "#A69E33", "#734002", "#595622", "#8C873F",  "#474B71")
colm = c("genome_a", "genome_b", "ani_true")
df_meta <- read.delim("all_pairs.tsv", comment.char = "#", header = FALSE, col.names=colm, colClasses = c("character", "character", "numeric")) |>
  mutate(
    pair = paste(pmin(genome_a, genome_b), pmax(genome_a, genome_b), sep = "|"),
    # ani_bin = cut(ani_true, breaks = c(55, 60, 65, 70, 75, 80, 85, 90, 95, 97.5, 99.9, 100), include.lowest = TRUE)
    ani_bin = cut(ani_true, breaks = c(55, 60, 80, 90, 95, 99, 100), include.lowest = TRUE)
  )
n_pairs <- nrow(df_meta)
df_meta <- df_meta %>% 
  mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
  mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
  select(!c(genome_a, genome_b)) %>%
  rename(genome_a = genome_x, genome_b = genome_y)
df_r <- rbind(
  vroom("results/ani-comparison/dashing2-estimates-all/distances/dashing2-v4.tsv.gz") %>%
    select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100),
  vroom("results/ani-comparison/mash-estimates-all/distances/mash-sensitive.tsv.gz") %>%
    select(method, genome_a, genome_b, ani_est = ani_pct, distance),
  vroom("results/ani-comparison/skani-estimates-all/distances/skani-slow.tsv.gz") %>%
    select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100),
  vroom("results/ani-comparison/fastani-estimates-all/distances/fastani-frag1000.tsv.gz") %>%
    select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100),
  vroom("results/ani-comparison/gdiff-estimates-all/distances/gdiff-k23-w23-h9-l333-n1000-frac50-chisq06635-p66.tsv.gz") %>%
    select(genome_a, genome_b, distance) %>% mutate(method = "gdiff", ani_est = (1-distance)*100)
) %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
  mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
  select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
df_r <- df_r %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
df <- merge(df_r, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)

df %>% filter(ani_true <= 99) %>%
  mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.25, 0.4), include.lowest = TRUE)) %>%
  ggplot() +
  # aes(x=ani_bin, fill=method, y=d_est/d_true, color=method) +
  aes(x=ani_bin, fill=method, y=((100-ani_est)-(100-ani_true))/(100-ani_true), color=method) +
  # aes(x=ani_bin, fill=method, y=((100-ani_est)-(100-ani_true))/(100-ani_true), color=method) +
  geom_boxplot(outliers = F, outlier.size = 0.6, position = position_dodge(0.8), color = "grey10") +
  # stat_summary(aes(group=method), geom="line") +
  # stat_summary(aes(group=method)) +
  # facet_wrap(~ani_bin, scale="free") +
  geom_hline(yintercept = 0, linetype=2) +
  # coord_cartesian(ylim = c(0.95, 1.15)) +
  # scale_y_log10() +
  theme_bw() + labs(x=TeX(r'(${D_{ANIb}$})'), y=TeX(r'(Bias)'), title="251,534 pairs from GORG", fill="Method") +
  scale_fill_manual(values = mc) +
  scale_color_manual(values = mc) + scale_y_continuous(labels=percent)
ggsave("./results/G-bias-lt99-anib.pdf", width = 7, height = 3.5)

mape_gt99 <- df %>%
  filter(ani_true > 99) %>%
  mutate(ani_bin = cut((100-ani_true)/100, c(0, 0.0025, 0.005, 0.01), include.lowest = T)) %>%
  group_by(ani_bin, method) %>%
  summarize(z=mean(abs((100-ani_true)-(100-ani_est))/(100-ani_true))*100, na.rm=T) %>%
  mutate(method=factor(method, levels=lvl)) %>%
  ggplot() +
  aes(fill=z, x=ani_bin, y=method) +
  geom_tile(show.legend = F) +
  geom_label(aes(label=round(z, 1), color=z>500), size = 2, show.legend = F, label.size = NA) +
  labs(fill=TeX(r'(MAPE)'), y="",
       title="MAPE",
       x=TeX(r'(${D_{ANIb}$})')) +
  # geom_abline(linetype="dashed") +
  theme_cowplot() +
  scale_fill_viridis_c(option = "B", direction = -1, values = c(0, 0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(values = c("black", "white")) +
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle=30, hjust = 1), legend.title = element_text(vjust=2))
mape_gt99
ggsave("./results/G-mape-gt99-anib.pdf", width = 1.85, height = 2.25)

mape_lte99 <- df %>%
  filter(ani_true <= 99) %>%
  filter(method!="skani" | ani_true > 80) %>%
  filter(method!="dashing2" | ani_true > 75) %>%
  mutate(ani_bin = cut((100-ani_true)/100, c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35), include.lowest = T)) %>%
  group_by(ani_bin, method) %>%
  summarize(z=mean(abs((100-ani_true)-(100-ani_est))/(100-ani_true))*100, na.rm=T) %>%
  complete(z) %>%
  mutate(method=factor(method, levels=lvl)) %>%
  ggplot() +
  aes(fill=z, x=ani_bin, y=method) +
  geom_tile(show.legend = F) +
  geom_label(aes(label=round(z, 1), color=z>10), size = 2, show.legend = F, label.size = NA) +
  labs(fill=TeX(r'(MAPE)'), y="",
       title="MAPE",
       x=TeX(r'(${D_{ANIb}$})')) +
  # geom_abline(linetype="dashed") +
  geom_label(data = data.frame(
                 z = c(NA, NA, NA, NA, NA),  
                 method=c("skani", "skani", "skani", "dashing2", "dashing2"), 
                 ani_bin=c("(0.2,0.25]", "(0.25,0.3]", "(0.3,0.35]", "(0.25,0.3]", "(0.3,0.35]")
               ), aes(y=method, x=ani_bin, label="x"), size = 5, fill = NA, color="black", label.size=NA, size=3) +
  theme_classic() +
  scale_fill_viridis_c(option = "B", direction = -1, values = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(values = c("black", "white")) +
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle=30, hjust = 1), legend.title = element_text(vjust=2))
mape_lte99
ggsave("./results/G-mape-lte99-anib.pdf", width = 5.85, height = 2.25)

df %>% filter(ani_true <= 99) %>%
  # mutate(ani_bin = cut(ani_true, c(60, 75, 90, 95, 99))) %>%
  mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.25, 0.4), include.lowest = TRUE)) %>%
  ggplot() +
  # aes(x=ani_bin, fill=method, y=d_est/d_true, color=method) +
  # aes(x=ani_bin, fill=method, y=(100-ani_est)/(100-ani_true), color=method) +
  aes(x=method, fill=method, y=((100-ani_est)-(100-ani_true))/(100-ani_true), color=method) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), outliers = F, trim = T, scale = "width", outlier.size = 0.2, color = "grey10") +
  # stat_summary(aes(group=method), geom="line") +
  stat_summary(color="gray15", size=0.25) +
  facet_wrap(~ani_bin, scale="fixed", nrow=1) +
  geom_hline(yintercept = 0, linetype="dashed") +
  # coord_cartesian(ylim = c(0.95, 1.15)) +
  # scale_y_log10() +
  theme_bw() + labs( y=TeX(r'(Bias)'), title="251,534 pairs from GORG") +
  scale_color_manual(values = mc) + scale_y_continuous(labels=percent) +
  scale_fill_manual(values = mc) + theme(axis.text.x = element_blank()) +
  labs(x="", fill="Method") + coord_cartesian(ylim = c(-0.45, 0.60)) 
ggsave("./results/G-bias-violin-lt99-anib.pdf", width = 7, height = 3.25)

ggscatter(
  df %>% filter(ani_true > 90),
  x = "ani_true", y = "ani_est", color = "method",
  add = "reg.line", alpha=0.05, conf.int = TRUE
  ) +
  geom_abline(linewidth=1, alpha=1, linetype="dashed") +
  stat_cor(aes(color = method), method = "pearson", show.legend = F) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc) +
  coord_cartesian(x=c(90, 100), y=c(90, 100)) +
  labs(x="ANIb", y=TeX(r'(${\hat{ANI}}$)'), color="Method", fill="Method")

df %>% filter(ani_true > 90) %>%
  ggplot() +
  aes(x=ani_true, y=ani_est, color=method, fill=method) +
  stat_cor(aes(label = ..r.label..), method = "spearman", show.legend = F) +
  geom_point(alpha=0.05) +
  stat_smooth(method = "lm", linewidth=1.5) +
  geom_abline(linewidth=1, alpha=1, linetype="dashed") +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc) +
  theme_bw() +
  labs(x="ANIb", y=TeX(r'(${\hat{ANI}}$)'), color="Method", fill="Method") +
  coord_cartesian(x=c(90, 100), y=c(90, 100))

df %>%
  filter((ani_true < 90) & (ani_true > 80)) %>%
  # filter((ani_true > 90)) %>%
  ggplot() +
  aes(x=1-ani_true/100, y=1-ani_est/100, color=method, fill = method) +
  stat_cor(aes(label = ..r.label..), label.y = c(0.195, 0.190, 0.185, 0.180, 0.175), position = position_dodge2(), method = "pearson", show.legend = F) +
  geom_point(alpha=0.025) +
  stat_smooth(method = "lm", linewidth=1.2) +
  geom_abline(linewidth=1, alpha=1, linetype="dashed") +
  geom_text(
    data = function(d) d %>% 
      group_by(method) %>% 
      summarize(
        # mape = mean(abs(((100-ani_true) - (100-ani_est)) / (100-ani_true))) * 100, 
        mape = mean(abs(((100-ani_true) - (100-ani_est)) / (100-ani_true))) * 100, 
        .groups = "drop"
      ) %>% 
      mutate(
        x = 0.125,
        y = 0.1985 - ((row_number())* 0.005)
      ),
    aes(x=x, y=y, label=paste0("MAPE = ", round(mape, 2), "%")),
    hjust=0, vjust=0, show.legend = F, 
  ) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc) +
  theme_cowplot(font_size = 11) +
  coord_cartesian(x=c(0.1, 0.2), y=c(0.1, 0.2)) +
  labs(y=TeX(r'(${\hat{D}$})'), x=TeX(r'(${D_{ANIb}$})'), color="Method", fill="Method")
ggsave("./results/G-pcorr-mape-lt90gt80-anib.pdf", width = 5, height = 3.75)

df %>%
  filter((ani_true > 90)) %>%
  ggplot() +
  aes(x=1-ani_true/100, y=1-ani_est/100, color=method, fill = method) +
  stat_cor(aes(label = ..r.label..), label.y = c(0.095, 0.090, 0.085, 0.080, 0.075), position = position_dodge2(), method = "pearson", show.legend = F) +
  geom_point(alpha=0.025) +
  stat_smooth(method = "lm", linewidth=1.2) +
  geom_abline(linewidth=1, alpha=1, linetype="dashed") +
  geom_text(
    data = function(d) d %>% 
      group_by(method) %>% 
      summarize(
        # mape = mean(abs(((100-ani_true) - (100-ani_est)) / (100-ani_true))) * 100, 
        mape = mean(abs(((100-ani_true) - (100-ani_est)) / (100-ani_true))) * 100, 
        .groups = "drop"
      ) %>% 
      mutate(
        x = 0.025,
        y = 0.0985 - ((row_number())* 0.005)
      ),
    aes(x=x, y=y, label=paste0("MAPE = ", round(mape, 2), "%")),
    hjust=0, vjust=0, show.legend = F
  ) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc) +
  theme_cowplot(font_size = 11) +
  coord_cartesian(x=c(0.0, 0.1), y=c(0.0, 0.1)) +
  labs(y=TeX(r'(${\hat{D}$})'), x=TeX(r'(${D_{ANIb}$})'), color="Method", fill="Method")
ggsave("./results/G-pcorr-mape-gt90-anib.pdf", width = 5, height = 3.75)

df %>% filter(ani_true <= 99) %>%
  filter(method!="skani" | ani_true > 80) %>%
  filter(method!="dashing2" | ani_true > 75) %>%
  mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4), include.lowest = TRUE)) %>%
  group_by(ani_bin, method) %>%
  summarise(mae=mean(abs((100-ani_est)-(100-ani_true))/(100-ani_true)), se = sd(abs((100-ani_est)-(100-ani_true))/(100-ani_true), na.rm = TRUE) / sqrt(n())) %>%
  ggplot() +
  aes(x=ani_bin, y=mae, color=method) +
  geom_line(aes(group=method)) +
  geom_point(aes(shape=method), size=3, alpha=0.85) +
  geom_errorbar(
    aes(group=method, ymin = mae - se, ymax = mae + se), width = 0.2
  ) +
  # geom_boxplot(aes(x=ani_bin)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_bw() +
  scale_y_continuous(labels=percent) +
  labs(x=TeX(r'($D_{ANIb})'), y=TeX(r'(MAPE)'), 
       #title="251,534 pairs from GORG", 
       color="Method", shape="Method") +
  scale_fill_manual(values = mc) +
  scale_color_manual(values = mc) +
  scale_shape_manual(values = c(16, 16, 17, 16, 16)) +
  theme(axis.text.x = element_text(angle=30, hjust=1))
ggsave("./results/G-mape-summary-lt99-anib.pdf", width = 4, height = 2.25)

merge(
  df %>% mutate(ani_bin=cut(ani_true, c(65, 70, 75, 80, 100))) %>% group_by(method, ani_bin) %>% summarise(c=n()),
  df_meta %>% mutate(ani_bin=cut(ani_true, c(65, 70, 75, 80, 100))) %>% group_by(ani_bin) %>% summarise(t=n()), all.y=T
) %>% complete(method, ani_bin, fill = list(c = 0, t = 1)) %>%
  ggplot() +
  aes(x=ani_bin, fill=method, y=c/t) +
  geom_col(color="gray20", position = position_dodge2(width = 0.8, padding = 0.15), na.rm = F) +
  # geom_line(aes(group=method, color=method)) +
  theme_bw() + labs(x="ANIb", y="Pairs with an estimate (%)", title="251,534 pairs from GORG") +
  scale_fill_manual(values = mc) +
  # scale_color_manual(values = mc) +
  scale_y_continuous(labels = percent)
ggsave("./results/G-completeness_percent-bar-anib.pdf", width = 4, height = 2.5)

cnpairs <- df %>%
  group_by(ani_true, method) %>% summarise(nx=n()) %>%
  group_by(method) %>%
  arrange(ani_true) %>%
  mutate(nt = cumsum(nx)) %>%
  # mutate(ani_true=round(1-ani_true/100, 3)) %>%
  mutate(ani_true=round(ani_true, 1)) %>%
  group_by(ani_true, method) %>%
  summarize(nt=mean(nt)) %>%
  ggplot() +
  aes(x=1-(1-ani_true/100), y=nt, color=method) +
  geom_line(aes(group=method, linetype = method), linewidth=2) +
  geom_hline(yintercept = 251534, color="darkgray", linetype="dashed") +
  # stat_smooth(aes(group=method), m) +
  theme_bw() +
  # labs(x=TeX(r'($D_{ANIb})'), y="# of pairs", title="251,534 pairs from GORG") +
  labs(x=TeX(r'($D_{ANIb})'), y="Cumulative num. of pairs") + # title="251,534 pairs from GORG"
  scale_color_manual(values = mc) +
  scale_x_continuous(labels = c(0.3, 0.2, 0.1, 0.0)) +
  scale_linetype_manual(values = c(1, 1, 2, 3, 1))#  +
  # scale_linetype_manual(values = c(1, 1, 1, 1, 1))#  +
# scale_y_continuous(labels = percent)
cnpairs
ggsave("./results/G-completeness_count-cumulative-anib.pdf", width = 4, height = 2.5)

df <- vroom("gene_coordinates/blastn-comparison/genes_mapped.tsv")
names(df)[1] <- "Genome A"
names(df)[2] <- "Genome B"

df_ani <- df_ani <- vroom("anib-groundtruth.csv")
df_ani <- rbind(df_ani %>% filter(`Genome A` < `Genome B`), df_ani %>% filter(`Genome A` >= `Genome B`) %>% rename(`Genome A`=`Genome B`, `Genome B`=`Genome A`))

blastnc <- merge(df, df_ani) %>% 
  filter(blastn_aln_len > 200) %>%
  ggplot() +
  aes(x=blastn_d, y=gdiff_d, color=1-`ANIgenome, %`/100) +
  geom_point(alpha=0.25) +
  stat_cor(aes(label = ..r.label..), method = "pearson", label.y = 0.29) +
  geom_abline(color="red") +
  scale_color_viridis_c(end = 1, option = "G") +
  # geom_text(
  #  data = function(d) d %>% 
  #    filter(is.finite(blastn_d), is.finite(gdiff_d)) %>%
  #    summarize(mae = mean(abs(((blastn_d) - (gdiff_d)))),  .groups = "drop") %>% 
  #    mutate(x = 0.0, y = 0.25),
  #  aes(x=x, y=y, label=paste0("MAE = ", round(mae, 2), "")), hjust=0, vjust=0, show.legend = F, inherit.aes = F
  # ) +
  theme_bw() +
  panel_border(color = "black", size = 0.5) +
  theme(legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.4, "cm")) +
  coord_cartesian(xlim=c(0, 0.3), ylim=c(0, 0.3)) +
  # labs(x="BLASTn", y=TeX(r'(gdiff $D$)'), color="ANI", title="Distances between annotated genes in marine microbes")
  labs(x="BLASTn", y="gdiff", color=TeX(r'($D_{ANIb})')) #, title="Distances between annotated genes")
blastnc
plot_grid(cnpairs + theme(legend.position = "none"),
          mape_gt99 + theme(),
          mape_lte99 + theme(axis.text.y = element_blank()) +
            labs(title="") + theme(axis.text.x.bottom = element_text(angle=30, hjust = 1), axis.title = element_text(size = 19)),
          blastnc,
          labels = c('', '', ''),
          label_size = 12,
          rel_widths = c(3.75, 2.5, 3.25, 4),
          nrow = 1)
ggsave("./results/gorg-ani-summary.pdf", width = 12, height = 3)

# dashing2 configuration comparison
if (T) {
df_dashing <- rbind(
  vroom("results/ani-comparison/dashing2-estimates-all/distances/dashing2-v1.tsv.gz") %>% mutate(method="v1"),
  vroom("results/ani-comparison/dashing2-estimates-all/distances/dashing2-v2.tsv.gz") %>% mutate(method="v2"),
  vroom("results/ani-comparison/dashing2-estimates-all/distances/dashing2-v3.tsv.gz") %>% mutate(method="v3"),
  vroom("results/ani-comparison/dashing2-estimates-all/distances/dashing2-v4.tsv.gz") %>% mutate(method="v4"),
  vroom("results/ani-comparison/dashing2-estimates-all/distances/dashing2-v5.tsv.gz") %>% mutate(method="v5")
) %>% select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100) %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
  mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
  select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
df_dashing <- df_dashing %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
df_dashing <- merge(df_dashing, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)
df_dashing %>% filter(ani_true <= 99) %>%
  mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4), include.lowest = TRUE)) %>%
  ggplot() +
  aes(x=ani_bin, color=method, y=abs((100-ani_est)-(100-ani_true))/(100-ani_true)) +
  stat_summary(aes(group=method), geom="line") +
  stat_summary() +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme_bw() + labs(x=TeX(r'(${D_{ANIb}$})'), y=TeX(r'(MAPE)'), title="251,534 pairs from GORG", color="Method") +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(labels=percent) +
  theme(axis.text.x = element_text(angle=30, hjust=1))
  ggsave("./results/G-mape-dashing2-config-anib.pdf", width = 5, height = 3.5)
}

# mash configuration comparison
if (T) {
  df_mash <- rbind(
    vroom("results/ani-comparison/mash-estimates-all/distances/mash-sensitive.tsv.gz") %>% mutate(method="sensitive"),
    vroom("results/ani-comparison/mash-estimates-all/distances/mash-long-k.tsv.gz") %>% mutate(method="long-k"),
    vroom("results/ani-comparison/mash-estimates-all/distances/mash-large-sketch.tsv.gz") %>% mutate(method="large-sketch"),
    vroom("results/ani-comparison/mash-estimates-all/distances/mash-default.tsv.gz") %>% mutate(method="default")
  ) %>% select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100) %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
    mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
    select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
  df_mash <- df_mash %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
  df_mash <- merge(df_mash, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)
  df_mash %>% filter(ani_true <= 99) %>%
    mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4), include.lowest = TRUE)) %>%
    ggplot() +
    aes(x=ani_bin, color=method, y=abs((100-ani_est)-(100-ani_true))/(100-ani_true)) +
    stat_summary(aes(group=method), geom="line") +
    stat_summary() +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_hline(yintercept = 0, linetype="dashed") +
    theme_bw() + labs(x=TeX(r'(${D_{ANIb}$})'), y=TeX(r'(MAPE)'), title="251,534 pairs from GORG", color="Method") +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(labels=percent) +
    theme(axis.text.x = element_text(angle=30, hjust=1))
  ggsave("./results/G-mape-mash-config-anib.pdf", width = 5, height = 3.5)
}

# skani configuration comparison
if (T) {
  df_skani <- rbind(
    vroom("results/ani-comparison/skani-estimates-all/distances/skani-slow.tsv.gz") %>% mutate(method="slow"),
    vroom("results/ani-comparison/skani-estimates-all/distances/skani-sensitive.tsv.gz") %>% mutate(method="sensitive"),
    vroom("results/ani-comparison/skani-estimates-all/distances/skani-fast.tsv.gz") %>% mutate(method="fast"),
    vroom("results/ani-comparison/skani-estimates-all/distances/skani-default.tsv.gz") %>% mutate(method="default")
  ) %>% select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100) %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
    mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
    select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
  df_skani <- df_skani %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
  df_skani <- merge(df_skani, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)
  df_skani %>% filter(ani_true <= 99) %>%
    mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4), include.lowest = TRUE)) %>%
    ggplot() +
    aes(x=ani_bin, color=method, y=abs((100-ani_est)-(100-ani_true))/(100-ani_true)) +
    stat_summary(aes(group=method), geom="line") +
    stat_summary() +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_hline(yintercept = 0, linetype="dashed") +
    theme_bw() + labs(x=TeX(r'(${D_{ANIb}$})'), y=TeX(r'(MAPE)'), title="251,534 pairs from GORG", color="Method") +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(labels=percent) +
    theme(axis.text.x = element_text(angle=30, hjust=1))
  ggsave("./results/G-mape-skani-config-anib.pdf", width = 5, height = 3.5)
}

# fastani configuration comparison
if (T) {
  df_fastani <- rbind(
    vroom("results/ani-comparison/fastani-estimates-all/distances/fastani-frag1000.tsv.gz") %>% mutate(method="--fragLen 1000 --minFraction 0.1"),
    vroom("results/ani-comparison/fastani-estimates-all/distances/fastani-frag3000.tsv.gz") %>% mutate(method="--fragLen 3000 --minFraction 0.1")
  ) %>% select(method, genome_a, genome_b, ani_est = ani_pct) %>% mutate(distance = 1-ani_est/100)  %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
    mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
    select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
  df_fastani <- df_fastani %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
  df_fastani <- merge(df_fastani, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)
  df_fastani %>% filter(ani_true <= 99) %>%
    mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4), include.lowest = TRUE)) %>%
    ggplot() +
    aes(x=ani_bin, color=method, y=abs((100-ani_est)-(100-ani_true))/(100-ani_true)) +
    stat_summary(aes(group=method), geom="line") +
    stat_summary() +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_hline(yintercept = 0, linetype="dashed") +
    theme_bw() + labs(x=TeX(r'(${D_{ANIb}$})'), y=TeX(r'(MAPE)'), title="251,534 pairs from GORG", color="Method") +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(labels=percent) +
    theme(axis.text.x = element_text(angle=30, hjust=1))
  ggsave("./results/G-mape-fastani-config-anib.pdf", width = 5, height = 3.5)
}

# gdiff configuration comparison
if (T) {
  df_gdiff <- rbind(
    vroom("results/ani-comparison/gdiff-estimates-all/distances/gdiff-k23-w23-h10-l500-n1000-frac50-chisq06635-p66.tsv.gz") %>% mutate(method="k23-w23-h10-l500-n1000-frac50"),
    vroom("results/ani-comparison/gdiff-estimates-all/distances/gdiff-k23-w23-h11-l500-n1000-frac50-chisq03841-p66.tsv.gz") %>% mutate(method="k23-w23-h11-l500-n1000-frac50"),
    vroom("results/ani-comparison/gdiff-estimates-all/distances/gdiff-k23-w23-h9-l333-n1000-frac50-chisq06635-p66.tsv.gz") %>% mutate(method="k23-w23-h9-l333-n1000-frac50"),
    vroom("results/ani-comparison/gdiff-estimates-all/distances/gdiff-k23-w23-h9-l500-n1000-frac33-chisq10828-p66.tsv.gz") %>% mutate(method="k23-w23-h9-l500-n1000-frac33"),
    vroom("results/ani-comparison/gdiff-estimates-all/distances/gdiff-k23-w23-h9-l500-n1000-frac50-chisq10828-p66.tsv.gz") %>% mutate(method="k23-w23-h9-l500-n1000-frac50")
    ) %>% select(method, genome_a, genome_b, distance) %>% mutate(ani_est = 100-distance*100) %>% mutate(genome_x = if_else(genome_a > genome_b, genome_b, genome_a)) %>%
    mutate(genome_y = if_else(genome_a > genome_b, genome_a, genome_b)) %>%
    select(!c(genome_a, genome_b)) %>% rename(genome_a = genome_x, genome_b = genome_y)
  df_gdiff <- df_gdiff %>% group_by(genome_a, genome_b, method) %>% summarise(ani_est=mean(ani_est, na.rm = T), distance=mean(distance, na.rm = T))
  df_gdiff <- merge(df_gdiff, df_meta) %>% mutate(d_true = 1 - ani_true/100) %>% rename(d_est = distance)
  df_gdiff %>% filter(ani_true <= 99) %>%
    mutate(ani_bin=cut(1-ani_true/100, breaks = c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4), include.lowest = TRUE)) %>%
    ggplot() +
    aes(x=ani_bin, color=method, y=abs((100-ani_est)-(100-ani_true))/(100-ani_true)) +
    stat_summary(aes(group=method), geom="line") +
    stat_summary() +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_hline(yintercept = 0, linetype="dashed") +
    theme_bw() + labs(x=TeX(r'(${D_{ANIb}$})'), y=TeX(r'(MAPE)'), title="251,534 pairs from GORG", color="Method") +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(labels=percent) +
    theme(axis.text.x = element_text(angle=30, hjust=1))
  ggsave("./results/G-mape-gdiff-config-anib.pdf", width = 5, height = 3.5)
}

