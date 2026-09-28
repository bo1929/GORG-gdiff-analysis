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

# detach("package:fitdistrplus", unload = TRUE)
# detach("package:MASS", unload = TRUE)
make_symmetric_max <- function(df) {
  df <- rbind(df %>% filter(genome_a > genome_b), df %>% filter(genome_b > genome_a) %>%rename(genome_a=genome_b, genome_b=genome_a)) %>%
    group_by(genome_a, genome_b) %>% summarize(ani_pct=max(ani_pct, na.rm=T))
  df <- df %>% filter(grepl("gnd", genome_a)) 
  return(df)
}
make_symmetric_mean <- function(df) {
  df <- rbind(df %>% filter(genome_a > genome_b), df %>% filter(genome_b > genome_a) %>%rename(genome_a=genome_b, genome_b=genome_a)) %>%
    group_by(genome_a, genome_b) %>% summarize(ani_pct=mean(ani_pct, na.rm=T))
  df <- df %>% filter(grepl("gnd", genome_a)) 
  return(df)
}

df_fastani <- vroom("results/distances/fastani-frag1000.tsv") %>% 
  select(genome_a, genome_b, ani_pct)
df_fastani <- make_symmetric_max(df_fastani) %>% mutate(method="fastANI")
df_mash <- vroom("results/distances/mash-k19-s10000.tsv") %>% 
  select(genome_a, genome_b, ani_pct) %>% mutate(method="mash") %>%
  filter(grepl("gnd", genome_a))
df_skani <- vroom("results/distances/skani-slow-min-af0.tsv") %>% 
  select(genome_a, genome_b, ani_pct) %>% mutate(method="skani") %>%
  filter(grepl("gnd", genome_a))
df_dashing2 <- vroom("results/distances/dashing2-v4.tsv") %>% 
  select(genome_a, genome_b, ani_pct) %>% mutate(method="dashing2") %>%
  filter(grepl("gnd", genome_a))
df_gdiff <- vroom("results/distances/gdiff-frac50-k23-w23-h9-l333-n1000-delta3-chisq06635-p66.tsv") %>% 
  mutate(ani_pct=100-100*distance) %>% 
  select(genome_a, genome_b, ani_pct) %>%
  mutate(method="gdiff")

pairs <- df_mash %>% select(genome_a, genome_b)
dfm <- vroom("metadata.tsv") %>% rename(genome_b=seed, genome_a=key)
dfm <- merge(pairs, dfm)
dfm <- dfm  %>% 
  mutate(missing_percent = if_else(is.na(missing_percent), 0, missing_percent)) %>%
  mutate(s = if_else(is.na(s), 0, s)) %>%
  mutate(ratevar=if_else(grepl("a22", level), "low", "high")) %>%
  mutate(true_ani_bin=cut(true_ani, c(80, 90, 95, 99, 100))) %>%
  mutate(missing_percent_bin=cut(missing_percent, c(0, 1, 5, 10, 15, 20), include.lowest = T)) %>%
  mutate(true_ani_bin=cut(true_ani, c(80, 90, 95, 99, 100))) %>%
  mutate(missing_percent_bin=cut(missing_percent, c(0, 1, 5, 10, 15, 20), include.lowest = T))

df <- merge(rbind(df_mash, df_skani, df_gdiff, df_fastani, df_dashing2), dfm)
df 

# df %>% filter(missing_percent == 0) %>%
#   filter(ratevar == "high") %>%
#   ggplot() +
#   aes(x=true_ani, y=ani_pct, color=method) +
#   geom_point() +
#   theme_bw() + labs(x="True ANI", y="Estimated")
mc = c("#809D6F", "#D0D55C", "#BD4030", "#C3A97E", "#769DA6", "#A69E33", "#734002", "#595622", "#8C873F",  "#474B71")

df %>%
  mutate(true_ani_bin=cut(1-true_ani/100, c(0, 0.01, 0.05, 0.1, 0.2))) %>%
  filter(ratevar == "high") %>%
  filter(true_ani != 100) %>%
  # filter (missing_percent <=15) %>%
  ggplot() +
  facet_wrap(~true_ani_bin, nrow=1) +
  # facet_grid(missing_percent_bin~true_ani_bin, scales = "free_y") +
  # aes(x=variation, y=(ani_pct)/(true_ani), fill=method) +
  # geom_hline(yintercept = 0, linetype = 14) +
  aes(x=missing_percent_bin, y=abs((100-ani_pct)-(100-true_ani))/(100-true_ani)) +
  # geom_boxplot(outliers = F) +
  stat_summary(aes(color=method)) +
  stat_summary(geom="line", aes(group=method, color=method), alpha=0.66, show.legend = FALSE, position = position_dodge2()) +
  stat_summary(aes(color=method), fun.data = mean_se, geom = "errorbar", width = 0.3) +
  # labs(y=TeX(r'(${\hat{D}}/{D}$)'), x="Missing portion (%)", color="Method") +
  labs(y=TeX(r'(Mean Absolute Error (%))'), x="Missing portion (%)", color="Method") +
  theme_bw() +
  scale_y_log10() +
  scale_y_continuous(labels = percent, trans="log", breaks = c(0.001, 0.01, 0.1, 1)) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc)
ggsave("./S-mape-lt100-missing_data-long.pdf", width = 11, height = 2)
# ggsave("./S-mape-lt100-missing_data.pdf", width = 8, height = 4)
# We should clarify in the label of dashing2 that it's containment.

df %>%
  mutate(true_ani_bin=cut(1-true_ani/100, c(0.01, 0.05, 0.1, 0.2))) %>%
  # filter(ratevar == "high") %>%
  filter(true_ani != 100) %>%
  # filter (missing_percent <=15) %>%
  filter (missing_percent == 0) %>%
  ggplot() +
  # facet_grid(missing_percent_bin~true_ani_bin, scales = "free_y") +
  # aes(x=variation, y=(ani_pct)/(true_ani), fill=method) +
  # geom_hline(yintercept = 0, linetype = 14) +
  aes(x=method, y=abs((100-ani_pct)-(100-true_ani))/(100-true_ani)) +
  # geom_boxplot(outliers = F) +
  geom_violin(aes(fill=method), draw_quantiles = c(0.25, 0.5, 0.75), trim = T) +
  # facet_wrap(~ratevar) +
  stat_summary(color="gray20") +
  # stat_summary(geom="line", aes(group=method, color=method), show.legend = FALSE, position = position_dodge2()) +
  # stat_summary(color="gray20", fun.data = mean_se, geom = "errorbar", width = 0.2) +
  # labs(y=TeX(r'(${\hat{D}}/{D}$)'), x="Missing portion (%)", color="Method") +
  labs(y=TeX(r'(Absolute Error (%))'), x="Methods", fill="Method") +
  theme_minimal_hgrid(font_size = 11) +
  # scale_y_log10() +
  scale_y_continuous(labels = percent, transform = "log", breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10)) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc) +
  coord_cartesian(ylim = c(0.0001, 1))
ggsave("./S-ape-lt100-violin.pdf", width = 6, height = 2)

df %>%
  mutate(true_ani_bin=cut(1-true_ani/100, c(0.01, 0.05, 0.1, 0.2))) %>%
  filter(ratevar == "high") %>%
  filter(true_ani <= 99) %>%
  ggplot() +
  # facet_grid(missing_percent_bin~true_ani_bin, scales = "free_y") +
  # aes(x=variation, y=(ani_pct)/(true_ani), fill=method) +
  geom_hline(yintercept = 0, linetype = 14) +
  aes(x=1-true_ani/100, y=((100-ani_pct)-(100-true_ani))/(100-true_ani), color=method) +
  geom_point(alpha=0.1, size=0.25) +
  stat_smooth(alpha=0.65, se=F) +
  labs(y=TeX(r'(Percentage Error)'), x="D", color="Method") +
  theme_minimal_grid(font_size=11) +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc)
ggsave("./S-pe-lt99-stat_smooth-high-wmissing.pdf", width = 5, height = 3)

df %>%
  mutate(true_ani_bin=cut(1-true_ani/100, c(0.01, 0.05, 0.1, 0.2))) %>%
  filter(ratevar == "low") %>%
  filter(true_ani <= 99) %>%
  ggplot() +
  # facet_grid(missing_percent_bin~true_ani_bin, scales = "free_y") +
  # aes(x=variation, y=(ani_pct)/(true_ani), fill=method) +
  geom_hline(yintercept = 0, linetype = 14) +
  aes(x=1-true_ani/100, y=((100-ani_pct)-(100-true_ani))/(100-true_ani), color=method) +
  geom_point(alpha=0.1, size=0.25) +
  stat_smooth(alpha=0.65, se=F) +
  labs(y=TeX(r'(Percentage Error)'), x="D", color="Method") +
  theme_minimal_grid(font_size=11) +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc)
ggsave("./S-pe-lt99-stat_smooth-rate-complete.pdf", width = 5, height = 3)

df %>% filter(true_ani <= 99) %>%
  filter(missing_percent == 0) %>%
  mutate(true_ani_bin = 1-round(true_ani)/100) %>%
  group_by(method, ratevar, true_ani_bin) %>%
  summarise(mean_error=mean(((100-ani_pct)-(100-true_ani))/(100-true_ani))) %>%
  pivot_wider(names_from = ratevar, values_from = mean_error) %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 14) +
  aes(color=method) +
  geom_segment(aes(x = true_ani_bin+as.numeric(as.factor(method))/1000, xend = true_ani_bin+as.numeric(as.factor(method))/1000, y = low, yend=high), arrow = grid::arrow(type="closed", length = unit(5, "pt"))) +
  labs(y=TeX(r'(Mean Percentage Error)'), x="D", color="Method") +
  theme_minimal_grid(font_size=11) +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc)
ggsave("./S-mpe-lt99-ratevar_change.pdf", width = 5, height = 3)

df %>% filter(true_ani <= 99) %>%
  filter(ratevar == "high") %>%
  # filter(missing_percent == 0 | missing_percent == 10) %>%
  mutate(wmissing = missing_percent > 0) %>%
  mutate(true_ani_bin = 1-round(true_ani)/100) %>%
  group_by(method, wmissing, true_ani_bin) %>%
  summarise(mean_error=mean(((100-ani_pct)-(100-true_ani))/(100-true_ani))) %>%
  pivot_wider(names_from = wmissing, values_from = mean_error) %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = 14) +
  aes(color=method) +
  geom_segment(aes(x = true_ani_bin+as.numeric(as.factor(method))/1000, xend = true_ani_bin+as.numeric(as.factor(method))/1000, y = `FALSE`, yend=`TRUE`), arrow = grid::arrow(type="closed", length = unit(5, "pt"))) +
  labs(y=TeX(r'(Mean Percentage Error)'), x="D", color="Method") +
  theme_minimal_grid(font_size=11) +
  scale_y_continuous(labels = percent) +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc)
ggsave("./S-mpe-lt99-missing_change.pdf", width = 5, height = 3)

df %>%
  filter(ratevar == "high") %>%
  mutate(true_ani_bin=cut(1-true_ani/100, c(0, 0.01, 0.05, 0.1, 0.2))) %>%
  filter(true_ani !=100) %>%
  filter (missing_percent == 0) %>%
  ggplot() +
  facet_wrap(~true_ani_bin, scales = "free_y") +
  # facet_grid(missing_percent_bin~true_ani_bin, scales = "free_y") +
  # aes(x=variation, y=(ani_pct)/(true_ani), fill=method) +
  geom_hline(yintercept = 0, linetype = 14) +
  aes(x=method, y=(true_ani-ani_pct)/100) +
  geom_boxplot(aes(fill=method), outliers = T) +
  # stat_summary(aes(color=method)) +
  # stat_summary(geom="line", aes(group=method, color=method), show.legend = FALSE, position = position_dodge2()) +
  # labs(y=TeX(r'(${\hat{D}}/{D}$)'), x="Missing portion (%)", color="Method") +
  labs(y=TeX(r'(${\hat{D}}-{D}$)'), fill="Method", x="Method") +
  theme_bw() +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc)
ggsave("./S-bias-lt100-high-complete.pdf", width = 8, height = 4.5)

df %>%
  filter(missing_percent == 0) %>%
  filter(true_ani != 100) %>%
  mutate(true_ani_bin=cut(1-true_ani/100, c(0, 0.01, 0.1, 0.2))) %>%
  ggplot() +
  facet_wrap(~true_ani_bin, scales = "fixed") +
  aes(x=reorder(ratevar, ani_pct/true_ani), y=abs((100-ani_pct)-(100-true_ani))/(100-true_ani), color=method) +
  geom_hline(yintercept = 0, linetype = 14) +
  stat_summary(geom="line", alpha=0.66, aes(group=method)) +
  # stat_summary(aes(shape = method)) +
  stat_summary(aes(shape=method)) +
  # geom_point(alpha=0.1, size=0.15) +
  # geom_violin(alpha=0.1, size=0.15) +
  # labs(y=TeX(r'(${\hat{D}}/{D}$)'), x="Rate variation") +
  labs(y=TeX(r'(Mean Absolute Error (%))'), x="Rate variation", color="Method", shape="Method") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) +
  scale_color_manual(values = mc) +
  scale_color_manual(values = mc) +
  scale_y_continuous(labels=percent) +
  scale_shape_manual(values = c(16, 16, 17, 16, 16))
ggsave("./S-mape-lt100-ratevar_change-alt.pdf", width = 6, height = 3)

df %>%
  filter(true_ani != 100) %>%
  filter(missing_percent == 0) %>%
  filter(ratevar == "high") %>%
  ggplot() +
  # facet_wrap(~ratevar) +
  aes(x=(100-true_ani)/100, y=(100-ani_pct)/100, color=method, fill=method) +
  geom_point(alpha=0.15) +
  stat_smooth() +
  labs(y=TeX(r'(${\hat{D}}$)'), x=TeX(r'(${{D}}$)'), color="Method", fill="Method") +
  geom_abline(linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc) +
  # coord_cartesian(xlim=c(0, 0.18), ylim=c(0, 0.18))+
  facet_zoom(x = (100-true_ani)/100 < 0.05, y = (100-ani_pct)/100 < 0.05, , zoom.size = 1, shrink = T) +
  geom_text(
    data = function(d) d %>% 
      group_by(method) %>% 
      summarize(
        mape = mean(abs(((100-true_ani) - (100-ani_pct)) / (100-true_ani))) * 100, 
        .groups = "drop"
      ) %>% 
      mutate(
        x = 0.06,
        y = 0.19 - ((row_number())* 0.01)
      ),
    aes(x=x, y=y, label=paste0("MAPE = ", round(mape, 2), "%")),
    hjust=0, vjust=0, show.legend = F, size = 3
  )
ggsave("./S-point-zoom-high.pdf", width = 7, height = 3.5)

df %>%
  filter(missing_percent == 0) %>% # filter(true_ani > 99) %>%
  ggplot() +
  facet_wrap(~ratevar) +
  aes(x=(100-true_ani)/100, y=(100-ani_pct)/100, color=method, shape=method) +
  geom_point(alpha=0.1) +
  stat_summary_bin() +
  labs(y=TeX(r'(${\hat{D}}$)'), x=TeX(r'(${{D}}$)')) +
  geom_abline(linetype="dashed") +
  theme_bw() +
  scale_color_manual(values = mc) +
  scale_fill_manual(values = mc) +
  facet_zoom(x = (100-true_ani)/100 < 0.05, y = (100-ani_pct)/100 < 0.05, , zoom.size = 1, shrink = T)

df %>%
  filter(true_ani != 100) %>%
  filter(missing_percent == 0) %>% filter(true_ani > 99) %>%
  mutate(true_ani_bin = cut((100-true_ani)/100, c(0, 0.0025, 0.005, 0.01), include.lowest = T)) %>%
  group_by(true_ani_bin, method, ratevar) %>%
  summarize(z=mean(abs((100-true_ani)-(100-ani_pct))/((100-true_ani)))*100) %>%
  ggplot() +
  facet_wrap(~ratevar) +
  aes(fill=z, y=true_ani_bin, x=method) +
  geom_tile() +
  geom_label(aes(label=round(z, 1), color=z>5), show.legend = F) +
  labs(fill=TeX(r'(MAPE)'), x="Method", y=TeX(r'(${D_{ANIb}$})')) +
  # geom_abline(linetype="dashed") +
  theme_cowplot(font_size = 10) +
  scale_fill_viridis_c(option = "B", direction = -1) +
  scale_color_manual(values = c("black", "white")) +
  theme(axis.text.x.bottom = element_text(angle=0.45))
ggsave("./S-mae-gt99.pdf", width = 8, height = 2.5)

mc = c("#809D6F", "#D0D55C", "#BD4030", "#C3A97E", "#769DA6", "#A69E33", "#734002", "#595622", "#8C873F",  "#474B71")

# Comparing different configurations of gdiff
if (F) {
df_gdiff <- merge(
  rbind(
    vroom("results/distances/gdiff-frac50-k23-w23-h9-l500-n1000-delta3-chisq10828-p66.tsv") %>%  mutate(method="gdiff: -h 9 -l 500 --frac 0.5"),
    vroom("results/distances/gdiff-frac50-k23-w23-h10-l500-n1000-delta3-chisq10828-p66.tsv") %>%  mutate(method="gdiff: -h 10 -l 500 --frac 0.5"),
    vroom("results/distances/gdiff-frac50-k23-w23-h11-l500-n1000-delta3-chisq10828-p66.tsv") %>%  mutate(method="gdiff: -h 11 -l 500 --frac 0.5"),
    vroom("results/distances/gdiff-frac33-k23-w23-h9-l500-n1000-delta3-chisq10828-p66.tsv") %>%  mutate(method="gdiff: -h 9 -l 500 --frac 0.33"),
    vroom("results/distances/gdiff-frac50-k23-w23-h9-l333-n1000-delta3-chisq06635-p66.tsv") %>%  mutate(method="gdiff (default): -h 9 -l 333 --frac 0.5  (*)")
  )%>% 
    mutate(ani_pct=100-100*distance) %>% 
    select(genome_a, genome_b, ani_pct, method),
  dfm
)
df_gdiff %>%
  mutate(true_ani_bin=cut(1-true_ani/100, c(0, 0.01, 0.1, 0.2))) %>%
  filter(ratevar == "high") %>% filter(true_ani != 100) %>%
  filter(missing_percent == 0) %>%
  ggplot() +
  geom_hline(yintercept = 0.00001, linetype = 14) +
  aes(x=true_ani_bin, y=abs(true_ani-ani_pct)/(100-true_ani), fill=method) +
  geom_boxplot(outliers = F) +
  # stat_summary(aes(color=method)) +
  # stat_summary(geom="line", aes(group=method, color=method), position = position_dodge2()) +
  labs(y=TeX(r'(Absolute error (%))'), x="Missing portion (%)", fill="Configuration") +
  theme_bw() +
  scale_y_continuous(transform = "log", labels=percent, breaks=c(0.0001, 0.001, 0.01, 0.1)) +
  scale_color_manual(values=c("#BD4030", "#9D4030", "#7D4030", "#5D4030", "#3D4030")) +
  scale_fill_manual(values=c("#BD4030", "#9D4030", "#7D4030", "#5D4030", "#3D4030")) +
  theme(legend.position = "bottom", legend.direction = "vertical")
ggsave("./S-gdiff-config-compare.pdf", width = 5, height = 4)
}

# Comparing different configurations of skani
if (T) {
  df_skani <- merge(
    rbind(
      vroom("results/distances/skani-slow-min-af0.tsv") %>%  mutate(method="skani: --slow --min-af 0 (*)"),
      vroom("results/distances/skani-c125-min-af0.tsv") %>%  mutate(method="skani: default"), 
      vroom("results/distances/skani-robust-min-af0-c70.tsv") %>%  mutate(method="skani: --robust -c 70")
    )%>% 
      select(genome_a, genome_b, ani_pct, method),
    dfm
  )
  df_skani %>%
    mutate(true_ani_bin=cut(1-true_ani/100, c(0, 0.01, 0.1, 0.2))) %>%
    filter(ratevar == "high") %>% filter(true_ani != 100) %>%
    filter(missing_percent == 0) %>%
    ggplot() +
    geom_hline(yintercept = 0.00001, linetype = 14) +
    aes(x=true_ani_bin, y=abs(true_ani-ani_pct)/(100-true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    # stat_summary(aes(color=method)) +
    # stat_summary(geom="line", aes(group=method, color=method), position = position_dodge2()) +
    labs(y=TeX(r'(Absolute error (%))'), x="Missing portion (%)", fill="Configuration") +
    theme_bw() +
    scale_y_continuous(transform = "log", labels=percent, breaks=c(0.0001, 0.001, 0.01, 0.1)) +
    scale_fill_manual(values=c("#668DCA", "#968DCA", "#769DA6")) +
    theme(legend.position = "bottom", legend.direction = "vertical")
ggsave("./S-skani-config-compare.pdf", width = 4, height = 4)
}

# Comparing different configurations of mash
if (T) {
  df_mash <- merge(
    rbind(
      vroom("results/distances/mash-k19-s10000.tsv") %>% mutate(method="mash: -k 19 -s 10000 (*)"), 
      vroom("results/distances/mash-k21-s1000.tsv") %>% mutate(method="mash: -k 21 -s 1000"), 
      vroom("results/distances/mash-k21-s10000.tsv") %>% mutate(method="mash: -k 21 -s 10000"), 
      vroom("results/distances/mash-k23-s10000.tsv") %>% mutate(method="mash: -k 23 -s 10000")
    )%>% 
      mutate(ani_pct=100-100*distance) %>% 
      select(genome_a, genome_b, ani_pct, method),
    dfm
  )
  df_mash %>%
    mutate(true_ani_bin=cut(1-true_ani/100, c(0, 0.01, 0.1, 0.2))) %>%
    filter(ratevar == "high") %>% filter(true_ani != 100) %>%
    filter(missing_percent == 0) %>%
    ggplot() +
    geom_hline(yintercept = 0.00001, linetype = 14) +
    aes(x=true_ani_bin, y=abs(true_ani-ani_pct)/(100-true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    # stat_summary(aes(color=method)) +
    # stat_summary(geom="line", aes(group=method, color=method), position = position_dodge2()) +
    labs(y=TeX(r'(Absolute error (%))'), x="Missing portion (%)", fill="Configuration") +
    theme_bw() +
    scale_y_continuous(transform = "log", labels=percent, breaks=c(0.0001, 0.001, 0.01, 0.1)) +
    scale_fill_manual(values=c("#C3A97E", "#C1B95E", "#C1991E", "#C17900")) +
    theme(legend.position = "bottom", legend.direction = "vertical")
  ggsave("./S-mash-config-compare.pdf", width = 5, height = 4)
}

# Comparing different configurations of fastani
if (T) {
  df_fastani <- merge(
    rbind(
      vroom("results/distances/fastani-frag1000.tsv") %>%  mutate(method="fastani: --fragLen 1000 --minFraction 0.1 (*)"), 
      vroom("results/distances/fastani-frag3000.tsv") %>%  mutate(method="fastani --fragLen 3000 --minFraction 0.1")
    )%>% 
      # mutate(ani_pct=100-100*distance) %>% 
      select(genome_a, genome_b, ani_pct, method),
    dfm
  )
  df_fastani %>%
    mutate(true_ani_bin=cut(1-true_ani/100, c(0, 0.01, 0.1, 0.2))) %>%
    filter(ratevar == "high") %>% filter(true_ani != 100) %>%
    filter(missing_percent == 0) %>%
    ggplot() +
    geom_hline(yintercept = 0.00001, linetype = 14) +
    aes(x=true_ani_bin, y=abs(true_ani-ani_pct)/(100-true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    # stat_summary(aes(color=method)) +
    # stat_summary(geom="line", aes(group=method, color=method), position = position_dodge2()) +
    labs(y=TeX(r'(Absolute error (%))'), x="Missing portion (%)", fill="Configuration") +
    theme_bw() +
    scale_y_continuous(transform = "log", labels=percent, breaks=c(0.0001, 0.001, 0.01, 0.1)) +
    scale_fill_manual(values=c("#EFB111", "#D0D55C")) +
    theme(legend.position = "bottom", legend.direction = "vertical")
  ggsave("./S-fastani-config-compare.pdf", width = 5, height = 4)
}

# Comparing different configurations of dashing2
if (T) {
  df_dashing2 <- merge(
    rbind(
      vroom("results/distances/dashing2-v4.tsv") %>%  mutate(method="dashing2: --symmetric-containment -k 23 -S 2048 (*)"), 
      vroom("results/distances/dashing2-default.tsv") %>%  mutate(method="dashing2: --mash-distance -k 31 -S 1024"), 
      vroom("results/distances/dashing2-containment.tsv") %>%  mutate(method="dashing2: --containment")
    )%>% 
      mutate(ani_pct=100-100*distance) %>% 
      select(genome_a, genome_b, ani_pct, method),
    dfm
  )
  df_dashing2%>%
    mutate(true_ani_bin=cut(1-true_ani/100, c(0, 0.01, 0.1, 0.2))) %>%
    filter(ratevar == "high") %>% filter(true_ani != 100) %>%
    filter(missing_percent == 0) %>%
    ggplot() +
    geom_hline(yintercept = 0.00001, linetype = 14) +
    aes(x=true_ani_bin, y=abs(true_ani-ani_pct)/(100-true_ani), fill=method) +
    geom_boxplot(outliers = F) +
    # stat_summary(aes(color=method)) +
    # stat_summary(geom="line", aes(group=method, color=method), position = position_dodge2()) +
    labs(y=TeX(r'(Absolute error (%))'), x="Missing portion (%)", fill="Configuration") +
    theme_bw() +
    scale_y_continuous(transform = "log", labels=percent, breaks=c(0.0001, 0.001, 0.01, 0.1)) +
    scale_fill_manual(values=c("#D0D55C", "#A09D6F", "#809D6F")) +
    theme(legend.position = "bottom", legend.direction = "vertical")
  ggsave("./S-dashing2-config-compare.pdf", width = 5, height = 4)
}
