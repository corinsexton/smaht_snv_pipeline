#!/usr/bin/env Rscript
# Prevalence of CrossTech / CrossCaller / CrossCore / CrossTissue per sample
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("~/smaht/smahtSNV_v2_core_specific/smaht_snv_pipeline/tests")

FINAL_DIR <- "../results_new_class/13_final"
OUT_FILE  <- "cross_tags.pdf"
BCFTOOLS  <- "/home/cos689/software/bin/bcftools"

vcfs <- list.files(FINAL_DIR, pattern = "\\.final\\.vcf\\.gz$", full.names = TRUE)
if (length(vcfs) == 0) stop("No .final.vcf.gz files found in ", FINAL_DIR)

TAGS <- c("CrossTech", "CrossCaller", "CrossCore", "CrossTissue")

dat <- lapply(vcfs, function(vcf) {
  sample_id <- sub("\\.final\\.vcf\\.gz$", "", basename(vcf))
  fmt <- paste(paste0("%", TAGS), collapse = " ")
  raw <- system2(BCFTOOLS,
                 c("query", "-f", paste0('"', fmt, '\\n"'), vcf),
                 stdout = TRUE)
  if (length(raw) == 0) return(NULL)
  df <- read.table(text = raw, header = FALSE, sep = " ",
                   col.names = TAGS, na.strings = ".")
  # Convert: 1 -> TRUE, NA -> FALSE
  for (tag in TAGS) df[[tag]] <- !is.na(df[[tag]])
  df$Sample <- sample_id
  df
})
dat <- bind_rows(dat)

# Summary: count and proportion of each tag per sample
tag_summary <- dat %>%
  group_by(Sample) %>%
  summarise(
    Total       = n(),
    across(all_of(TAGS), list(n = sum), .names = "{.col}_n"),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_n"),
                list(pct = ~ . / Total * 100),
                .names = "{.col}ct")) %>%
  rename_with(~ sub("_n_npct$", "_pct", .), ends_with("_npct"))

# Long form for plotting
long_n <- tag_summary %>%
  select(Sample, Total, ends_with("_n")) %>%
  pivot_longer(ends_with("_n"),
               names_to = "Tag", values_to = "Count") %>%
  mutate(Tag = sub("_n$", "", Tag),
         Pct = Count / Total * 100)

# Order samples by total variant count
sample_order <- tag_summary %>% arrange(Total) %>% pull(Sample)
long_n$Sample <- factor(long_n$Sample, levels = sample_order)
long_n$Tag    <- factor(long_n$Tag, levels = TAGS)

tag_colors <- c(
  CrossTech    = "#2196F3",
  CrossCaller  = "#4CAF50",
  CrossCore    = "#FF9800",
  CrossTissue  = "#9C27B0"
)

# --- Plot 1: raw counts, grouped bars ---
p_counts <- ggplot(long_n, aes(x = Sample, y = Count, fill = Tag)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = tag_colors) +
  labs(title = "Cross* tag counts per sample",
       x = NULL, y = "Variant count", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# --- Plot 2: percentage, grouped bars ---
p_pct <- ggplot(long_n, aes(x = Sample, y = Pct, fill = Tag)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = tag_colors) +
  labs(title = "Cross* tag prevalence per sample (% of all variants)",
       x = NULL, y = "% of variants", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# --- Plot 3: heatmap — % of variants with each tag ---
p_heat <- ggplot(long_n, aes(x = Tag, y = Sample, fill = Pct)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%", Pct)), size = 3) +
  scale_fill_gradient(low = "white", high = "#1565C0",
                      name = "% variants") +
  labs(title = "Cross* tag prevalence heatmap",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# --- Plot 4: upset-style co-occurrence — stacked bar of flag combinations ---
combo_dat <- dat %>%
  mutate(combo = paste0(
    ifelse(CrossTech,    "Tech+",    ""),
    ifelse(CrossCaller,  "Caller+",  ""),
    ifelse(CrossCore,    "Core+",    ""),
    ifelse(CrossTissue,  "Tissue+",  "")
  )) %>%
  mutate(combo = ifelse(combo == "", "None", sub("\\+$", "", combo))) %>%
  count(Sample, combo) %>%
  group_by(Sample) %>%
  mutate(Pct = n / sum(n) * 100) %>%
  ungroup()

combo_dat$Sample <- factor(combo_dat$Sample, levels = sample_order)

# Order combos by overall frequency within each CrossTissue group,
# then place CrossTissue-absent first, CrossTissue-present second
combo_freq <- combo_dat %>%
  group_by(combo) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total))

has_tissue  <- grepl("Tissue", combo_freq$combo)
combo_order <- c(combo_freq$combo[!has_tissue], combo_freq$combo[has_tissue])
combo_dat$combo <- factor(combo_dat$combo, levels = combo_order)

n_no  <- sum(!has_tissue)
n_yes <- sum(has_tissue)

no_tissue_colors <- if (n_no  > 0) colorRampPalette(c("#b3d9f7", "#08467a"))(n_no)  else character(0)
tissue_colors    <- if (n_yes > 0) colorRampPalette(c("#e8c6f0", "#5b0070"))(n_yes) else character(0)

combo_colors <- setNames(c(no_tissue_colors, tissue_colors), combo_order)

p_combo <- ggplot(combo_dat, aes(x = Sample, y = Pct, fill = combo)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = combo_colors, name = "Flag combination") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Cross* flag combination proportions per sample",
       x = NULL, y = "% of variants") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

pdf(OUT_FILE, width = 14, height = 6)
print(p_counts)
ggsave("cross_tags_counts.png", plot = p_counts, height = 4, width = 7, dpi = 300)
print(p_pct)
ggsave("cross_tags_pct.png", plot = p_pct, height = 4, width = 7, dpi = 300)
print(p_combo)
ggsave("cross_tags_combo.png", plot = p_combo, height = 4, width = 7, dpi = 300)
dev.off()

pdf(sub("\\.pdf$", "_heatmap.pdf", OUT_FILE), width = 8, height = 7)
print(p_heat)
ggsave("cross_tags_heatmap.png", plot = p_heat, height = 5, width = 5, dpi = 300)
dev.off()

message("Wrote: ", OUT_FILE)
message("Wrote: ", sub("\\.pdf$", "_heatmap.pdf", OUT_FILE))
