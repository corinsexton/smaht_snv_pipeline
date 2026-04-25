#!/usr/bin/env Rscript
# EvidenceScore distribution across all 13_final VCFs
library(ggplot2)
library(dplyr)

setwd("~/smaht/smahtSNV_v2_core_specific/smaht_snv_pipeline/tests")

FINAL_DIR <- "../results_new_class/13_final"
OUT_FILE  <- "evidence_score_distribution.pdf"
BCFTOOLS  <- "/home/cos689/software/bin/bcftools"

vcfs <- list.files(FINAL_DIR, pattern = "\\.final\\.vcf\\.gz$", full.names = TRUE)
if (length(vcfs) == 0) stop("No .final.vcf.gz files found in ", FINAL_DIR)

dat <- lapply(vcfs, function(vcf) {
  sample_id <- sub("\\.final\\.vcf\\.gz$", "", basename(vcf))
  raw <- system2(BCFTOOLS,
                 c("query", "-f", '"%FILTER %EvidenceScore\\n"', vcf),
                 stdout = TRUE)
  if (length(raw) == 0) return(NULL)
  df <- read.table(text = raw, header = FALSE, sep = " ",
                   col.names = c("Filter", "EvidenceScore"))
  df$Sample <- sample_id
  df
})
dat <- bind_rows(dat)
dat$EvidenceScore <- factor(dat$EvidenceScore, levels = 0:4)

# Order samples by total PASS count (ascending) for a meaningful layout
pass_order <- dat %>%
  filter(Filter == "PASS") %>%
  count(Sample) %>%
  arrange(n) %>%
  pull(Sample)
dat$Sample <- factor(dat$Sample, levels = pass_order)

score_colors <- c(
  "0" = "#d62728",   # LowEvidence (EvidenceScore = 0)
  "1" = "#aec7e8",
  "2" = "#1f77b4",
  "3" = "#0d4a8a",
  "4" = "#041e42"
)

# --- Plot 1: stacked bar (counts) ---
p_counts <- ggplot(dat, aes(x = Sample, fill = EvidenceScore)) +
  geom_bar(width = 0.75) +
  scale_fill_manual(values = score_colors,
                    name = "EvidenceScore",
                    labels = c("0 (LowEvidence)", "1", "2", "3", "4 (all flags)")) +
  labs(title = "EvidenceScore distribution per sample",
       x = NULL, y = "Variant count") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# --- Plot 2: stacked bar (proportions) ---
p_prop <- ggplot(dat, aes(x = Sample, fill = EvidenceScore)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_manual(values = score_colors,
                    name = "EvidenceScore",
                    labels = c("0 (LowEvidence)", "1", "2", "3", "4 (all flags)")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "EvidenceScore distribution per sample (proportions)",
       x = NULL, y = "Proportion") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# --- Plot 3: faceted bar per sample ---
p_facet <- ggplot(dat, aes(x = EvidenceScore, fill = EvidenceScore)) +
  geom_bar(width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = score_colors) +
  facet_wrap(~ Sample, scales = "free_y", ncol = 4) +
  labs(title = "EvidenceScore distribution — per sample",
       x = "EvidenceScore", y = "Variant count") +
  theme_bw(base_size = 8) +
  theme(strip.text = element_text(size = 6))

print(p_counts)
print(p_prop)
ggsave("proportion_ev_score.png",plot = p_prop,height = 4, width = 5, dpi = 300)

dev.off()

pdf(sub("\\.pdf$", "_facet.pdf", OUT_FILE), width = 14, height = 10)
print(p_facet)
ggsave("facet_ev_score.png",plot = p_facet,height = 6, width = 7, dpi = 300)

dev.off()

message("Wrote: ", OUT_FILE)
message("Wrote: ", sub("\\.pdf$", "_facet.pdf", OUT_FILE))
