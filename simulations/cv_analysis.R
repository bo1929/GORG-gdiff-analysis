library(vroom)
library(ggplot2)
library(dplyr)
# install.packages("ggpmisc")
library(ggpmisc)

LR_UB_THR   <- 3.841
MIN_PORTION <- 0.66

select_sample <- function(d, lr_ub, lr_ub_thr = LR_UB_THR, min_portion = MIN_PORTION) {
  kepp <- is.finite(d)
  d <- d[kepp]; lr_ub <- lr_ub[kepp]
  if (!length(d)) return(numeric(0))
  pass <- !is.na(lr_ub) & lr_ub > lr_ub_thr
  if (sum(pass) / length(d) > min_portion) d[pass] else d
}

pair_ani <- function(d, lr_ub) {
  sel <- select_sample(d, lr_ub)
  dplyr::tibble(
    ani_pct      = if (length(sel)) 100 - 100 * mean(sel) else NA_real_,
    stdev_pident = if (length(sel) > 1) 100 * sd(sel) else NA_real_
  )
}

read_samples <- function(path) {
  vroom(path, delim = "\t", show_col_types = FALSE) %>%
    filter(grepl("a5", genome_a), grepl("_p", genome_a)) %>%
    mutate(genome_a = gsub("_a5", "", genome_a),
           lr_ub = suppressWarnings(as.numeric(lr_ub)))
}

# ---------------- Load metadata ----------------
meta <- vroom("metadata.tsv", delim = "\t", show_col_types = FALSE)

# a5 keys, pruned (_p), with true ANI / missing% (like eval-simulations.R)
dfm <- meta %>%
  filter(grepl("a5", key), grepl("_p", key)) %>%
  mutate(genome_b = seed,
         genome_a_key = gsub("_a5", "", key)) %>%
  mutate(missing_percent = if_else(is.na(missing_percent), 0, missing_percent)) %>%
  dplyr::select(genome_a_key, genome_b, true_ani, missing_percent)

# ---------------- BLAST (results-ORFvANI/results-combined.csv) ----------------
# genome_b = base (_p0_s0), genome_a = variant (_pNN_sNN). Strip the "_c<x>"
# suffix / "_a5" marker so genome_a matches metadata's genome_a_key.
blast <- vroom("results-ORFvANI/results-combined.csv", delim = ",", show_col_types = FALSE) %>%
  filter(grepl("_p0_s0", genome_b)) %>%
  mutate(
    genome_b = gsub("_gnd00000_.*", "", genome_b),
    genome_a = gsub("_c.*", "", genome_a)
  ) %>%
  mutate(genome_a = gsub("_a5", "", genome_a)) %>%
  transmute(genome_a, genome_b, ani_pct = mean_pident, stdev_pident, method = "BLAST")

# ---------------- gdiff samples (lr_ub present -> filtered) ----------------
gdiff <- read_samples("results/samples/gdiff-frac50-k23-w23-h9-l500-n1000-delta3.tsv") %>%
  group_by(genome_a, genome_b) %>%
  summarise(pair_ani(d, lr_ub), .groups = "drop") %>%
  mutate(method = "gdiff")

# ---------------- fastANI samples (no lr_ub -> unfiltered) ----------------
fastani <- read_samples("results/samples/fastani-frag1000.tsv") %>%
  group_by(genome_a, genome_b) %>%
  summarise(pair_ani(d, lr_ub), .groups = "drop") %>%
  mutate(method = "fastANI")

# ---------------- Combine with metadata (CV) ----------------
cv <- bind_rows(
  gdiff,
  fastani,
  blast
) %>%
  inner_join(dfm, by = c(genome_a = "genome_a_key", genome_b = "genome_b")) %>%
  filter(true_ani != 100) %>%
  mutate(cv_pct = stdev_pident / (100 - ani_pct) * 100)

# ---------------- Figure ----------------
gg <- cv %>% filter(missing_percent == 0) %>%
  mutate(true_ani_bin = cut(true_ani, c(80, 99, 100), include.lowest = TRUE)) %>%
  ggplot() +
  aes(x = true_ani, y = cv_pct, color = method) +
  geom_point(alpha = 0.4) +
  facet_wrap(~true_ani_bin, scale = "free") +
  geom_hline(yintercept = 1 / sqrt(5) * 100, linetype = "dashed") +
  stat_smooth(method = "lm") +
  stat_poly_eq(aes(label = after_stat(eq.label)), formula = y ~ x, parse = TRUE) +
  theme_bw() +
  labs(x = "ANI", y = "Coefficient of variance (%)") +
  scale_color_manual(values = c("#474B71", "#D0D55C", "#9D4030"))

gg
