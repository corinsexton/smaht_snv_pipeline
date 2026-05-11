library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(forcats)

setwd("~/smaht/smahtSNV_v2_core_specific/smaht_snv_pipeline/tests")
# ── Parse z file ──────────────────────────────────────────────────────────────
z_path  <- file.path(dirname(getwd()), "tests/z")   # adjust if needed
if (!file.exists(z_path)) z_path <- "z"

lines <- readLines(z_path)

header_idx   <- grep("^=== ", lines)
summary_idx  <- grep("V2 total:", lines)

extract_num <- function(line, label) {
  m <- regmatches(line, regexpr(paste0(label, ":\\s*([0-9]+)"), line))
  if (length(m) == 0) return(NA_integer_)
  as.integer(sub(paste0(label, ":\\s*"), "", m))
}

records <- lapply(seq_along(header_idx), function(i) {
  sample_id <- sub("^=== (.+) ===$", "\\1", lines[header_idx[i]])
  s <- summary_idx[summary_idx > header_idx[i]][1]
  if (is.na(s)) return(NULL)
  ln <- lines[s]
  data.frame(
    sample    = sample_id,
    v2_total  = extract_num(ln, "V2 total"),
    v1_total  = extract_num(ln, "V1 total"),
    shared    = extract_num(ln, "Shared"),
    unique_v2 = extract_num(ln, "Unique V2"),
    unique_v1 = extract_num(ln, "Unique V1")
  )
})

df <- bind_rows(records) |>
  mutate(
    donor  = sub("-.*", "", sample),
    tissue = sub("^[^-]+-", "", sample)
  )

# ── Load tissue metadata ───────────────────────────────────────────────────────
meta <- read_tsv("~/smaht/metadata/coverages_tissues.tsv", show_col_types = FALSE) |>
  select(tissue_layer, tissue_desc, tissue) |>
  distinct() |>
  filter(!is.na(tissue), tissue != "")

# join on tissue code
df <- df |>
  left_join(meta, by = "tissue") |>
  mutate(
    tissue_desc  = if_else(is.na(tissue_desc), tissue, tissue_desc),
    tissue_layer = if_else(is.na(tissue_layer), "Unknown", tissue_layer)
  )

# ── Sort: by tissue_layer then tissue_desc ────────────────────────────────────
layer_order <- c("Clin", "Endo", "Meso", "Ecto", "Germ", "Unknown")

df <- df |>
  mutate(
    tissue_layer = factor(tissue_layer, levels = layer_order),
    tissue_desc  = factor(tissue_desc,
                          levels = df |>
                            arrange(tissue_layer, tissue_desc) |>
                            pull(tissue_desc) |>
                            unique())
  ) |>
  arrange(tissue_layer, tissue_desc, donor)

# ── Core counts per donor+tissue ──────────────────────────────────────────────
cores_raw <- read_tsv(
  "~/smaht/analysis/snv_25_donors_paper/aliquots_cores/ill_pb.tsv",
  col_names = FALSE, show_col_types = FALSE
)
core_counts <- cores_raw |>
  select(donor = X1, tissue = X2, core = X3) |>
  distinct(donor, tissue, core) |>
  count(donor, tissue, name = "n_cores")

df <- df |> left_join(core_counts, by = c("donor", "tissue")) |>
  mutate(n_cores = replace_na(n_cores, 0L))

# label: donor + tissue_desc + core count for x axis
df <- df |>
  mutate(label = paste0(donor, "\n", tissue_desc, " (", n_cores, ")"))

df$label <- factor(df$label, levels = unique(df$label))

# ── Reshape to long for ggplot ─────────────────────────────────────────────────
df_long <- df |>
  select(label, tissue_layer, tissue_desc, unique_v1, shared, unique_v2) |>
  pivot_longer(c(unique_v1, shared, unique_v2),
               names_to = "category", values_to = "count") |>
  mutate(category = factor(category,
                           levels = c("unique_v1", "shared", "unique_v2"),
                           labels = c("Unique V1", "Shared", "Unique V2")))

# ── Plot ───────────────────────────────────────────────────────────────────────
pal <- c("Unique V1" = "#d62728", "Shared" = "#7f7f7f", "Unique V2" = "#1f77b4")

p <- ggplot(df_long, aes(x = label, y = count, fill = category)) +
  geom_col(width = 0.75, position = "dodge") +
  facet_grid(. ~ tissue_layer, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "V1 vs V2 variant calls per tissue sample",
    x     = NULL,
    y     = "Variant count"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "#e8e8e8"),
    # strip.text       = element_text(face = "bold"),
    legend.position  = "top",
    panel.grid.major.x = element_blank()
  )

plot(p)
# ggsave("v1_v2_barplot.pdf", p, width = max(12, nrow(df) * 0.6), height = 6)
message("Saved: v1_v2_barplot.pdf")
