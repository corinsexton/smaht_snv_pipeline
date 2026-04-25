#!/usr/bin/env Rscript
# Cross* tag interaction and filter-field diagnostics
library(ggplot2) 
library(dplyr)
library(tidyr)

setwd("~/smaht/smahtSNV_v2_core_specific/smaht_snv_pipeline/tests")

FINAL_DIR <- "../results_new_class/13_final"
OUT_FILE  <- "cross_tag_interactions.pdf"
BCFTOOLS  <- "/home/cos689/software/bin/bcftools"

TAGS <- c("CrossTech", "CrossCaller", "CrossCore", "CrossTissue")

vcfs <- list.files(FINAL_DIR, pattern = "\\.final\\.vcf\\.gz$", full.names = TRUE)
if (length(vcfs) == 0) stop("No .final.vcf.gz files found in ", FINAL_DIR)

# Pull FILTER, EvidenceScore, and all 4 Cross* flags per variant per sample
dat <- lapply(vcfs, function(vcf) {
  sample_id <- sub("\\.final\\.vcf\\.gz$", "", basename(vcf))
  fmt <- paste(c("%FILTER", "%EvidenceScore", paste0("%", TAGS)), collapse = " ")
  raw <- system2(BCFTOOLS,
                 c("query", "-f", paste0('"', fmt, '\\n"'), vcf),
                 stdout = TRUE)
  if (length(raw) == 0) return(NULL)
  df <- read.table(text = raw, header = FALSE, sep = " ",
                   col.names = c("Filter", "EvidenceScore", TAGS),
                   na.strings = ".")
  for (tag in TAGS) df[[tag]] <- !is.na(df[[tag]])
  df$EvidenceScore <- as.integer(df$EvidenceScore)
  df$Sample <- sample_id
  df
})
dat <- bind_rows(dat)

# Sample order: ascending total variant count (consistent across plots)
sample_order <- dat %>% count(Sample) %>% arrange(n) %>% pull(Sample)
dat$Sample <- factor(dat$Sample, levels = sample_order)

tag_colors <- c(
  CrossTech   = "#2196F3",
  CrossCaller = "#4CAF50",
  CrossCore   = "#FF9800",
  CrossTissue = "#9C27B0"
)

###############################################################################
# Plot 1 — Pairwise conditional co-occurrence heatmap
#   Cell (row=A, col=B): P(B present | A present), across all variants/samples
###############################################################################
cond_mat <- outer(TAGS, TAGS, FUN = Vectorize(function(a, b) {
  rows_with_a <- dat[[a]]
  mean(dat[[b]][rows_with_a], na.rm = TRUE)
}))
rownames(cond_mat) <- TAGS
colnames(cond_mat) <- TAGS

cond_df <- as.data.frame(as.table(cond_mat)) %>%
  rename(TagA = Var1, TagB = Var2, Prob = Freq) %>%
  mutate(Label = sprintf("%.2f", Prob),
         TagA  = factor(TagA, levels = TAGS),
         TagB  = factor(TagB, levels = TAGS))

p_cond <- ggplot(cond_df, aes(x = TagB, y = TagA, fill = Prob)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#1565C0",
                      limits = c(0, 1), name = "P(col | row)") +
  labs(title = "Pairwise conditional co-occurrence",
       subtitle = "P(column tag present | row tag present), all samples pooled",
       x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

###############################################################################
# Plot 2 — Sole-rescuer breakdown (EvidenceScore == 1 variants)
#   Which single tag is solely responsible for PASS status, per sample?
###############################################################################
sole <- dat %>%
  filter(EvidenceScore == 1) %>%
  mutate(SoleTag = case_when(
    CrossTech    ~ "CrossTech",
    CrossCaller  ~ "CrossCaller",
    CrossCore    ~ "CrossCore",
    CrossTissue  ~ "CrossTissue",
    TRUE         ~ "Unknown"
  )) %>%
  count(Sample, SoleTag) %>%
  group_by(Sample) %>%
  mutate(Pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(SoleTag = factor(SoleTag, levels = TAGS))

p_sole <- ggplot(sole, aes(x = Sample, y = Pct, fill = SoleTag)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = tag_colors, name = "Sole rescuing tag") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Sole-rescuer breakdown (EvidenceScore = 1 PASS variants)",
       subtitle = "Which single tag saves each variant from LowEvidence?",
       x = NULL, y = "% of score-1 variants") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

###############################################################################
# Plot 3 — CrossTech stratification: EvidenceScore distribution by tech support
#   CrossTech+ (TIER1 proxy) vs CrossTech- (TIER2), stacked bars
###############################################################################
score_colors <- c(
  "0" = "#d62728",
  "1" = "#aec7e8",
  "2" = "#1f77b4",
  "3" = "#0d4a8a",
  "4" = "#041e42"
)

strat <- dat %>%
  mutate(TechGroup    = ifelse(CrossTech, "CrossTech+ (TIER1)", "CrossTech- (TIER2)"),
         EvidenceScore = factor(EvidenceScore, levels = 0:4))

p_strat_count <- ggplot(strat, aes(x = Sample, fill = EvidenceScore)) +
  geom_bar(width = 0.75) +
  scale_fill_manual(values = score_colors,
                    labels = c("0 (LowEvidence)", "1", "2", "3", "4"),
                    name = "EvidenceScore") +
  facet_wrap(~ TechGroup, ncol = 1, scales = "free_y") +
  labs(title = "EvidenceScore distribution stratified by CrossTech",
       x = NULL, y = "Variant count") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        strip.text = element_text(face = "bold"))

###############################################################################
# Plot 4 — Cumulative yield curve
#   For each min EvidenceScore threshold (0–4), % of variants that pass
#   One line per sample — shows sensitivity/stringency tradeoff
###############################################################################
thresholds <- 0:4

yield <- lapply(levels(dat$Sample), function(s) {
  sub_dat <- dat %>% filter(Sample == s)
  total   <- nrow(sub_dat)
  if (total == 0) return(NULL)
  data.frame(
    Sample    = s,
    Threshold = thresholds,
    Pct       = sapply(thresholds, function(t) sum(sub_dat$EvidenceScore >= t) / total * 100)
  )
})
yield <- bind_rows(yield)
yield$Sample <- factor(yield$Sample, levels = sample_order)

p_yield <- ggplot(yield, aes(x = Threshold, y = Pct, color = Sample, group = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = thresholds,
                     labels = paste0("≥", thresholds)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0, 100)) +
  labs(title = "Cumulative variant yield by EvidenceScore threshold",
       subtitle = "% of all variants with EvidenceScore ≥ threshold",
       x = "Minimum EvidenceScore", y = "% of variants retained",
       color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8))

###############################################################################
# Plot 5 — Tag co-occurrence rate per sample (line plot)
#   For each tag pair, what fraction of samples have both tags co-occurring
#   at rate > 0? Shows cross-sample consistency of co-occurrence.
###############################################################################
tag_pairs <- combn(TAGS, 2, simplify = FALSE)

cooccur_sample <- lapply(tag_pairs, function(pair) {
  a <- pair[1]; b <- pair[2]
  dat %>%
    group_by(Sample) %>%
    summarise(
      n_both  = sum(.data[[a]] & .data[[b]]),
      n_a     = sum(.data[[a]]),
      n_b     = sum(.data[[b]]),
      n_total = n(),
      .groups = "drop"
    ) %>%
    mutate(
      PairLabel = paste0(sub("Cross", "", a), " & ", sub("Cross", "", b)),
      Pct_of_total = n_both / n_total * 100
    )
}) %>% bind_rows()

cooccur_sample$Sample    <- factor(cooccur_sample$Sample,    levels = sample_order)
cooccur_sample$PairLabel <- factor(cooccur_sample$PairLabel)

pair_colors <- RColorBrewer::brewer.pal(length(tag_pairs), "Dark2")

p_cooccur <- ggplot(cooccur_sample, aes(x = Sample, y = Pct_of_total,
                                         color = PairLabel, group = PairLabel)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = pair_colors, name = "Tag pair") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Tag co-occurrence rate per sample",
       subtitle = "% of all variants where both tags in pair are present",
       x = NULL, y = "% variants with both tags") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

###############################################################################
# Plot 6 — Tag marginal prevalence with CI across samples
#   Box plot / dot plot of per-sample tag prevalence, one box per tag
###############################################################################
tag_prev_per_sample <- dat %>%
  group_by(Sample) %>%
  summarise(across(all_of(TAGS), ~ mean(.) * 100), .groups = "drop") %>%
  pivot_longer(all_of(TAGS), names_to = "Tag", values_to = "Pct") %>%
  mutate(Tag = factor(Tag, levels = TAGS))

p_box <- ggplot(tag_prev_per_sample, aes(x = Tag, y = Pct, fill = Tag)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = Tag), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_manual(values  = tag_colors, guide = "none") +
  scale_color_manual(values = tag_colors, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Cross* tag prevalence distribution across samples",
       subtitle = "Each point = one sample; box = median ± IQR",
       x = NULL, y = "% of variants with tag") +
  theme_bw(base_size = 12)

###############################################################################
# Write PDFs
###############################################################################
pdf(OUT_FILE, width = 12, height = 7)
print(p_cond)
ggsave("cross_interactions_cond.png",    plot = p_cond,        height = 5, width = 5,  dpi = 300)
print(p_sole)
ggsave("cross_interactions_sole.png",    plot = p_sole,        height = 4, width = 7,  dpi = 300)
print(p_strat_count)
ggsave("cross_interactions_strat.png",   plot = p_strat_count, height = 6, width = 7,  dpi = 300)
print(p_yield)
ggsave("cross_interactions_yield.png",   plot = p_yield,       height = 4, width = 7,  dpi = 300)
print(p_cooccur)
ggsave("cross_interactions_cooccur.png", plot = p_cooccur,     height = 4, width = 7,  dpi = 300)
print(p_box)
ggsave("cross_interactions_box.png",     plot = p_box,         height = 4, width = 5,  dpi = 300)
dev.off()

message("Wrote: ", OUT_FILE)
