library(data.table)
library(ggplot2)

data_in <- path.expand("../../data/variant_qual")
image_out <- path.expand("../../images/variant_qual")
if (!dir.exists(data_in)) stop("Data directory not found: ", data_in)
dir.create(image_out, recursive = TRUE, showWarnings = FALSE)
sources <- c("MCL_22", "MCL_586", "MCL_836")
keys <- c("CHROM", "POS", "REF", "ALT")

read_qc <- function(kind) {
  rbindlist(lapply(sources, function(source) {
    file <- file.path(data_in, sprintf("%s_quality_%s.tsv", kind, source))
    if (!file.exists(file)) stop("Missing output: ", file)
    x <- fread(file, na.strings = c(".", "NA"))
    x[, SOURCE := source]
    x
  }), use.names = TRUE, fill = TRUE)
}

sites <- read_qc("site")
carriers <- read_qc("carrier")

cat("Input records by source:\n")
print(sites[, .(
  n_sites = .N,
  n_KIR = sum(IS_KIR == 1L, na.rm = TRUE)
), by = SOURCE])

# Analyse only the three post_ppi MCL sources.
mcl_sites <- copy(sites)
mcl_carriers <- copy(carriers)

if (!nrow(mcl_sites)) stop("No post_ppi variant records found.")

# Audit variant records repeated across MCL files before deduplication.
site_duplicates <- mcl_sites[, .(
  n_records = .N,
  source_files = paste(unique(SOURCE), collapse = ";"),
  n_QUAL = uniqueN(QUAL, na.rm = TRUE),
  n_MQ = uniqueN(MQ, na.rm = TRUE),
  n_QD = uniqueN(QD, na.rm = TRUE),
  n_VQSLOD = uniqueN(VQSLOD, na.rm = TRUE),
  n_FS = uniqueN(FS, na.rm = TRUE),
  n_SOR = uniqueN(SOR, na.rm = TRUE),
  n_DP = uniqueN(DP, na.rm = TRUE),
  n_AC = uniqueN(AC, na.rm = TRUE),
  n_AF = uniqueN(AF, na.rm = TRUE),
  n_AN = uniqueN(AN, na.rm = TRUE),
  n_FILTER = uniqueN(FILTER, na.rm = TRUE),
  n_IS_KIR = uniqueN(IS_KIR, na.rm = TRUE)
), by = keys][n_records > 1L]

fwrite(site_duplicates, file.path(data_in, "mcl_repeated_sites.tsv"), sep = "\t")

checks <- setdiff(names(site_duplicates), c(keys, "n_records", "source_files"))
if (any(vapply(site_duplicates[, ..checks], function(x) any(x > 1L), logical(1)))) {
  stop("Conflicting quality or annotation for overlapping MCL sites. ",
       "Inspect mcl_repeated_sites.tsv before combining sources.")
}

# Preserve the union of KIR annotations and source provenance.
annotations <- mcl_sites[, .(
  IS_KIR = as.integer(any(IS_KIR == 1L)),
  KIR_GENES = {
    genes <- unique(unlist(strsplit(
      KIR_GENES[!is.na(KIR_GENES)], ";", fixed = TRUE
    )))
    if (length(genes)) paste(genes, collapse = ";") else NA_character_
  },
  SOURCE_FILES = paste(unique(SOURCE), collapse = ";")
), by = keys]

mcl_sites <- unique(mcl_sites, by = keys)
mcl_sites[, c("SOURCE", "IS_KIR", "KIR_GENES") := NULL]
mcl_sites <- merge(mcl_sites, annotations, by = keys, sort = FALSE)

# Audit duplicate sample/variant carrier observations.
call_keys <- c(keys, "SAMPLE")
call_duplicates <- mcl_carriers[, .(
  n_records = .N,
  n_GT = uniqueN(GT, na.rm = TRUE),
  n_DP = uniqueN(DP, na.rm = TRUE),
  n_GQ = uniqueN(GQ, na.rm = TRUE),
  n_AD_REF = uniqueN(AD_REF, na.rm = TRUE),
  n_AD_ALT = uniqueN(AD_ALT, na.rm = TRUE)
), by = call_keys][n_records > 1L]

fwrite(call_duplicates, file.path(data_in, "mcl_repeated_calls.tsv"), sep = "\t")

if (nrow(call_duplicates[
  n_GT > 1L | n_DP > 1L | n_GQ > 1L |
  n_AD_REF > 1L | n_AD_ALT > 1L
])) {
  stop("Conflicting sample-level calls across MCL files. ",
       "Inspect mcl_repeated_calls.tsv before combining sources.")
}

mcl_carriers <- unique(mcl_carriers, by = call_keys)
mcl_carriers[, SOURCE := NULL]
mcl_carriers <- merge(
  mcl_carriers,
  mcl_sites[, c(keys, "IS_KIR"), with = FALSE],
  by = keys, all.x = TRUE, sort = FALSE
)

if (anyNA(mcl_carriers$IS_KIR))
  stop("Some carrier calls have no matching MCL site annotation.")

if (mcl_sites[IS_KIR == 1L, .N] == 0L)
  stop("No KIR variants were annotated in any of the three MCL VCFs.")

mcl_sites[, region := fifelse(IS_KIR == 1L, "KIR", "Other post_ppi")]
mcl_carriers[, region := fifelse(IS_KIR == 1L, "KIR", "Other post_ppi")]

cat("Unique KIR sites:", mcl_sites[IS_KIR == 1L, .N], "\n")
print(mcl_sites[IS_KIR == 1L, .(
  CHROM, POS, REF, ALT, KIR_GENES, SOURCE_FILES
)])

fwrite(mcl_sites, file.path(data_in, "mcl_sites_deduplicated.tsv"), sep = "\t")

# Assign each distinct KIR site a plotting group. Sites annotated to several
# KIR genes are grouped as "Multiple KIR genes" rather than arbitrarily assigned.
qc_sites <- copy(mcl_sites)
qc_sites[, KIR_GENE := NA_character_]
qc_sites[IS_KIR == 1L, KIR_GENE := KIR_GENES]

if (qc_sites[IS_KIR == 1L, anyNA(KIR_GENE)]) {
  stop("Some KIR sites lack a KIR gene annotation.")
}

qc_sites[
  IS_KIR == 1L & grepl(";", KIR_GENE, fixed = TRUE),
  KIR_GENE := "Multiple KIR genes"
]

kir_genes <- sort(unique(qc_sites[IS_KIR == 1L, KIR_GENE]))
kir_colours <- setNames(
  grDevices::hcl.colors(length(kir_genes), palette = "Dark 3"),
  kir_genes
)

# Figure 1. Plot each chromosome 19 site from the MCL files as a point.
# No separately filtered chromosome 19 VCF is included.
is_chr19 <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE) == "19"

chr19_plot_sites <- unique(
  qc_sites[is_chr19(CHROM), .(
    CHROM, POS, REF, ALT, MQ, QD, VQSLOD, IS_KIR, KIR_GENE
  )],
  by = keys
)
if (!nrow(chr19_plot_sites)) stop("No chromosome 19 sites in MCL extracts.")
chr19_plot_sites[, Mb := POS / 1e6]

positional_metrics <- c("MQ", "QD", "VQSLOD")
positional_labels <- c(
  MQ = "Mapping quality (MQ)",
  QD = "Quality by depth (QD)",
  VQSLOD = "Variant quality log-odds (VQSLOD)"
)
positional_long <- melt(
  chr19_plot_sites,
  id.vars = c(keys, "Mb", "IS_KIR", "KIR_GENE"),
  measure.vars = positional_metrics,
  variable.name = "metric", value.name = "value"
)[is.finite(value)]
positional_long[, metric := factor(
  metric, levels = positional_metrics,
  labels = positional_labels[positional_metrics]
)]

p1 <- ggplot() +
  geom_point(
    data = positional_long[IS_KIR == 0L],
    aes(Mb, value),
    colour = "grey72", alpha = 0.45, size = 1.15
  ) +
  geom_point(
    data = positional_long[IS_KIR == 1L],
    aes(Mb, value, colour = KIR_GENE),
    alpha = 0.95, size = 2.2
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = kir_colours, name = "KIR gene") +
  labs(
    title = "Chromosome 19: individual variant quality",
    subtitle = "Grey: non-KIR post_ppi sites; coloured: KIR sites, all from the three MCL VCFs",
    x = "Chromosome 19 position (Mb)", y = "Observed quality metric"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(image_out, "kir_chr19_quality.png"), p1,
       width = 12, height = 9, dpi = 300)

kir_chr19_positions <- chr19_plot_sites[IS_KIR == 1L, POS]
if (length(kir_chr19_positions)) {
  zoom_min <- min(kir_chr19_positions)
  zoom_max <- max(kir_chr19_positions)
  zoom_padding <- max(250000, 0.25 * (zoom_max - zoom_min))
  p1_zoom <- p1 + coord_cartesian(
    xlim = c(max(0, zoom_min - zoom_padding), zoom_max + zoom_padding) / 1e6
  ) + labs(title = "Chromosome 19: quality near KIR genes")
  ggsave(file.path(image_out, "kir_chr19_quality_zoom.png"), p1_zoom,
         width = 12, height = 9, dpi = 300)
}

# Figure 2. Within-post_ppi site-quality comparison.
site_long <- melt(
  mcl_sites, id.vars = "region",
  measure.vars = c("MQ", "QD", "FS", "SOR", "VQSLOD"),
  variable.name = "metric", value.name = "value"
)[is.finite(value)]
site_long[, metric := factor(
  metric,
  levels = c("MQ", "QD", "FS", "SOR", "VQSLOD"),
  labels = c(
    "Mapping quality (MQ)",
    "Quality by depth (QD)",
    "Fisher strand bias (FS)",
    "Strand odds ratio (SOR)",
    "Variant quality log-odds (VQSLOD)"
  )
)]

p2 <- ggplot(site_long, aes(value, colour = region)) +
  stat_ecdf(geom = "step", linewidth = 0.7) +
  facet_wrap(~metric, scales = "free_x", ncol = 2) +
  labs(title = "KIR versus other post_ppi variants",
       x = "Quality metric", y = "Cumulative proportion", colour = NULL) +
  theme_bw(base_size = 11)

ggsave(file.path(image_out, "kir_site_quality_ecdf.png"), p2,
       width = 11, height = 8, dpi = 300)

# Figure 3. Within-post_ppi heterozygous carrier-quality comparison.
het <- mcl_carriers[GT %in% c("0/1", "1/0", "0|1", "1|0")]
het[, AB_deviation := abs(AB - 0.5)]

carrier_long <- melt(
  het, id.vars = "region",
  measure.vars = c("DP", "GQ", "AB_deviation"),
  variable.name = "metric", value.name = "value"
)[is.finite(value)]
carrier_long[, metric := factor(
  metric,
  levels = c("DP", "GQ", "AB_deviation"),
  labels = c(
    "Read depth (DP)",
    "Genotype quality (GQ)",
    "Allele balance deviation (|AB - 0.5|)"
  )
)]

if (nrow(carrier_long)) {
  p3 <- ggplot(carrier_long, aes(value, colour = region)) +
    stat_ecdf(geom = "step", linewidth = 0.7) +
    facet_wrap(~metric, scales = "free_x", ncol = 2) +
    labs(title = "KIR heterozygous carrier quality",
         x = "Quality metric", y = "Cumulative proportion", colour = NULL) +
    theme_bw(base_size = 11)
  
  ggsave(file.path(image_out, "kir_carrier_quality_ecdf.png"), p3,
         width = 11, height = 7, dpi = 300)
} else {
  message("No finite heterozygous carrier metrics to plot.")
}

site_summary <- mcl_sites[, .(
  n_sites = .N,
  pass_fraction = mean(FILTER == "PASS", na.rm = TRUE),
  median_MQ = median(MQ, na.rm = TRUE),
  median_QD = median(QD, na.rm = TRUE),
  median_VQSLOD = median(VQSLOD, na.rm = TRUE)
), by = region]

carrier_summary <- het[, .(
  n_heterozygous_calls = .N,
  median_DP = median(DP, na.rm = TRUE),
  median_GQ = median(GQ, na.rm = TRUE),
  median_AB_deviation = median(AB_deviation, na.rm = TRUE),
  fraction_AB_30_70 = mean(AB >= 0.3 & AB <= 0.7, na.rm = TRUE)
), by = region]

fwrite(site_summary, file.path(data_in, "kir_site_summary.tsv"), sep = "\t")
fwrite(carrier_summary, file.path(data_in, "kir_carrier_summary.tsv"), sep = "\t")
print(site_summary)
print(carrier_summary)

# Figures 4-8. Individual variants highlighted by KIR gene.
rank_metrics <- function(x) {
  x <- x[is.finite(value)]
  x[, percentile := 100 * (frank(value, ties.method = "average") - 0.5) / .N,
    by = metric]
  x
}

plot_ranked_metrics <- function(x, title, subtitle) {
  ggplot() +
    geom_vline(xintercept = c(5, 95), colour = "grey65",
               linetype = "dotted", linewidth = 0.35) +
    geom_point(
      data = x[IS_KIR == 0L],
      aes(percentile, value),
      colour = "grey75", alpha = 0.65, size = 1.2
    ) +
    geom_point(
      data = x[IS_KIR == 1L],
      aes(percentile, value, colour = KIR_GENE),
      alpha = 0.95, size = 2.3
    ) +
    facet_wrap(~metric, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = kir_colours, name = "KIR gene") +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    labs(
      title = title, subtitle = subtitle,
      x = "Percentile among analysed variants (per metric)",
      y = "Observed quality metric"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))
}

plot_histogram_with_kir_marks <- function(x, title, subtitle, bins = 35) {
  ggplot() +
    geom_histogram(
      data = x[IS_KIR == 0L],
      aes(x = value),
      bins = bins,
      fill = "grey80",
      colour = "white",
      linewidth = 0.25
    ) +
    geom_rug(
      data = x[IS_KIR == 1L],
      aes(x = value, colour = KIR_GENE),
      sides = "b",
      alpha = 0.95,
      linewidth = 0.8,
      length = grid::unit(0.08, "npc")
    ) +
    facet_wrap(~metric, scales = "free_x", ncol = 3) +
    scale_colour_manual(values = kir_colours, name = "KIR gene") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Observed quality metric",
      y = "Count of non-KIR variants"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))
}

# Figure 4. Every post_ppi site across six QC fields.
site_metrics <- c("MQ", "QD", "FS", "SOR", "VQSLOD", "QUAL")
site_metric_labels <- c(
  MQ = "Mapping quality (MQ)",
  QD = "Quality by depth (QD)",
  FS = "Fisher strand bias (FS)",
  SOR = "Strand odds ratio (SOR)",
  VQSLOD = "Variant quality log-odds (VQSLOD)",
  QUAL = "Variant quality score (QUAL)"
)
qc_long <- melt(
  qc_sites,
  id.vars = c(keys, "IS_KIR", "KIR_GENES", "KIR_GENE"),
  measure.vars = site_metrics,
  variable.name = "metric", value.name = "value"
)
qc_long <- rank_metrics(qc_long)
qc_long[, metric := factor(metric, levels = site_metrics, labels = site_metric_labels[site_metrics])]

p4 <- plot_ranked_metrics(
  qc_long,
  title = "Site quality across all post_ppi variants",
  subtitle = "Grey: all non-KIR sites; coloured: individual KIR genes. Dotted lines mark the 5th and 95th percentiles."
)

ggsave(file.path(image_out, "kir_all_site_quality_ranked.png"), p4,
       width = 14, height = 8.5, dpi = 300)

fwrite(
  qc_long[IS_KIR == 1L, .(
    CHROM, POS, REF, ALT, KIR_GENES, metric, value, percentile
  )],
  file.path(data_in, "kir_site_quality_percentiles.tsv"),
  sep = "\t"
)

# Figure 5. Two-dimensional view to reveal clustering and outliers.
qc_space <- qc_sites[is.finite(QD) & is.finite(VQSLOD)]
p5 <- ggplot() +
  geom_point(
    data = qc_space[IS_KIR == 0L],
    aes(QD, VQSLOD), colour = "grey75", size = 1.5, alpha = 0.60
  ) +
  geom_point(
    data = qc_space[IS_KIR == 1L],
    aes(QD, VQSLOD, colour = KIR_GENE), size = 2.7, alpha = 0.95
  ) +
  scale_colour_manual(values = kir_colours, name = "KIR gene") +
  labs(
    title = "Variant quality space: QD versus VQSLOD",
    subtitle = "One point per distinct post_ppi SNV; all non-KIR sites shown in grey",
    x = "Quality by depth (QD)", y = "Variant quality log-odds (VQSLOD)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(image_out, "kir_quality_space.png"), p5,
       width = 11, height = 8, dpi = 300)

# Figure 6. Median heterozygote-call quality per variant, not per sample.
het_sites <- het[, .(
  DP = median(DP, na.rm = TRUE),
  GQ = median(GQ, na.rm = TRUE),
  AB_deviation = median(AB_deviation, na.rm = TRUE)
), by = keys]

het_sites <- merge(
  het_sites,
  qc_sites[, c(keys, "IS_KIR", "KIR_GENE"), with = FALSE],
  by = keys, all.x = TRUE, sort = FALSE
)
if (anyNA(het_sites$IS_KIR)) stop("Unmatched heterozygous variant sites.")

het_metrics <- c("DP", "GQ", "AB_deviation")
het_metric_labels <- c(
  DP = "Read depth (DP)",
  GQ = "Genotype quality (GQ)",
  AB_deviation = "Allele balance deviation (|AB - 0.5|)"
)
het_long <- melt(
  het_sites,
  id.vars = c(keys, "IS_KIR", "KIR_GENE"),
  measure.vars = het_metrics,
  variable.name = "metric", value.name = "value"
)
het_long <- rank_metrics(het_long)
het_long[, metric := factor(metric, levels = het_metrics, labels = het_metric_labels[het_metrics])]

if (nrow(het_long)) {
  p6 <- plot_ranked_metrics(
    het_long,
    title = "Heterozygous carrier quality by variant",
    subtitle = "Per-variant medians. Lower allele-balance deviation is closer to the expected 0.5."
  )
  ggsave(file.path(image_out, "kir_carrier_quality_ranked.png"), p6,
         width = 14, height = 5.5, dpi = 300)
  
  # Figure 8. Histograms of per-variant carrier quality with KIR rug marks.
  p8 <- plot_histogram_with_kir_marks(
    het_long,
    title = "Distribution of heterozygous carrier quality by variant",
    subtitle = "Grey histograms show non-KIR variants. Coloured rug marks show individual KIR variants."
  )
  ggsave(file.path(image_out, "kir_carrier_quality_histograms.png"), p8,
         width = 14, height = 5.8, dpi = 300)
}

# Figure 7. Histograms of site-level quality with KIR rug marks.
p7 <- plot_histogram_with_kir_marks(
  qc_long,
  title = "Distribution of site quality across post_ppi variants",
  subtitle = "Grey histograms show non-KIR variants. Coloured rug marks show individual KIR variants."
)

ggsave(file.path(image_out, "kir_site_quality_histograms.png"), p7,
       width = 14, height = 8.5, dpi = 300)

# Figure 9 and statistical report. One distinct SNV is the analysis unit for
# site metrics; one per-variant median heterozygote call is the unit for carrier
# metrics. No carrier calls are treated as independent variants.
set.seed(20261016)

quality_metrics <- rbindlist(list(
  data.table(
    type = "Site",
    code = c("MQ", "QD", "FS", "SOR", "VQSLOD", "QUAL"),
    label = c(
      "Mapping quality (MQ)", "Quality by depth (QD)",
      "Fisher strand bias (FS)", "Strand odds ratio (SOR)",
      "Variant quality log-odds (VQSLOD)", "Variant quality (QUAL)"
    ),
    concerning = c("low", "low", "high", "high", "low", "low")
  ),
  data.table(
    type = "Carrier median",
    code = c("DP", "GQ", "AB_deviation"),
    label = c("Read depth (DP)", "Genotype quality (GQ)",
              "Allele balance deviation (|AB - 0.5|)"),
    concerning = c("low", "low", "high")
  )
))

# Two-sided, exploratory rank-sum tests. Cliff's delta describes stochastic
# separation (positive = higher values in KIR), independently of significance.
# Bootstrapped confidence intervals quantify uncertainty in the median difference
# (KIR minus non-KIR). These are site-resampling intervals, not read-level CIs.
compare_one_metric <- function(kir, other, boot_reps = 1500L) {
  kir <- kir[is.finite(kir)]
  other <- other[is.finite(other)]
  nk <- length(kir)
  no <- length(other)
  if (!nk || !no) {
    return(data.table(
      n_KIR = nk, n_other = no,
      median_KIR = NA_real_, median_other = NA_real_,
      median_difference = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      cliffs_delta = NA_real_, p_value = NA_real_
    ))
  }
  observed_diff <- median(kir) - median(other)
  if (nk >= 2L && no >= 2L) {
    boots <- replicate(
      boot_reps,
      median(sample(kir, nk, replace = TRUE)) -
        median(sample(other, no, replace = TRUE))
    )
    bounds <- as.numeric(quantile(boots, c(0.025, 0.975), names = FALSE))
  } else {
    bounds <- c(NA_real_, NA_real_)
  }
  if (nk >= 2L && no >= 2L &&
      !(length(unique(c(kir, other))) == 1L)) {
    test <- suppressWarnings(stats::wilcox.test(
      kir, other, alternative = "two.sided", exact = FALSE
    ))
    u <- as.numeric(test$statistic)
    effect <- 2 * u / (nk * no) - 1
    pvalue <- test$p.value
  } else if (length(unique(c(kir, other))) == 1L) {
    effect <- 0
    pvalue <- 1
  } else {
    effect <- NA_real_
    pvalue <- NA_real_
  }
  data.table(
    n_KIR = nk, n_other = no,
    median_KIR = median(kir), median_other = median(other),
    median_difference = observed_diff,
    ci_low = bounds[1], ci_high = bounds[2],
    cliffs_delta = effect, p_value = pvalue
  )
}

site_test_data <- copy(qc_sites)
carrier_test_data <- copy(het_sites)

stat_rows <- lapply(seq_len(nrow(quality_metrics)), function(i) {
  info <- quality_metrics[i]
  x <- if (info$type == "Site") site_test_data else carrier_test_data
  result <- compare_one_metric(
    x[IS_KIR == 1L][[info$code]],
    x[IS_KIR == 0L][[info$code]]
  )
  cbind(info, result)
})
qc_statistics <- rbindlist(stat_rows)
qc_statistics[, p_FDR := p.adjust(p_value, method = "BH")]
qc_statistics[, direction := fifelse(
  is.na(median_difference), "undetermined",
  fifelse(
    median_difference == 0, "same median",
    fifelse(
      (concerning == "low" & median_difference < 0) |
        (concerning == "high" & median_difference > 0),
      "KIR shifted in concerning direction",
      "KIR shifted in favourable direction"
    )
  )
)]
setcolorder(qc_statistics, c(
  "type", "code", "label", "concerning", "n_KIR", "n_other",
  "median_KIR", "median_other", "median_difference", "ci_low",
  "ci_high", "cliffs_delta", "p_value", "p_FDR", "direction"
))
fwrite(qc_statistics, file.path(data_in, "kir_qc_statistical_comparison.tsv"), sep = "\t")

# Count KIR observations lying in the theoretically concerning 5% tail
# of the empirical non-KIR distribution. Strict comparisons avoid marking
# every tied observation at the percentile boundary as an outlier.
metric_long_data <- function(x, info) {
  long <- melt(
    x,
    id.vars = c(keys, "IS_KIR", "KIR_GENE"),
    measure.vars = info$code,
    variable.name = "code", value.name = "value"
  )[is.finite(value)]
  long[, type := info$type[1L]]
  merge(long, info[, .(type, code, label, concerning)],
        by = c("type", "code"), sort = FALSE)
}

all_qc_values <- rbindlist(list(
  metric_long_data(site_test_data, quality_metrics[type == "Site"]),
  metric_long_data(carrier_test_data, quality_metrics[type == "Carrier median"])
), use.names = TRUE)

tail_thresholds <- all_qc_values[IS_KIR == 0L, .(
  tail_cutoff = if (concerning[1L] == "low")
    as.numeric(quantile(value, 0.05, names = FALSE))
  else as.numeric(quantile(value, 0.95, names = FALSE))
), by = .(type, code, concerning)]

all_qc_values <- merge(
  all_qc_values, tail_thresholds,
  by = c("type", "code", "concerning"), sort = FALSE
)
all_qc_values[, concerning_tail := ifelse(
  concerning == "low", value < tail_cutoff, value > tail_cutoff
)]
all_qc_values[, region := fifelse(IS_KIR == 1L, "KIR", "Other post_ppi")]

tail_summary <- all_qc_values[, .(
  n = .N,
  n_concerning_tail = sum(concerning_tail),
  fraction_concerning_tail = mean(concerning_tail),
  cutoff = tail_cutoff[1L]
), by = .(type, code, label, concerning, region)]

fwrite(tail_summary, file.path(data_in, "kir_qc_tail_summary.tsv"), sep = "\t")
fwrite(
  all_qc_values[IS_KIR == 1L & concerning_tail == TRUE, .(
    type, CHROM, POS, REF, ALT, KIR_GENE,
    code, label, value, concerning, tail_cutoff
  )],
  file.path(data_in, "kir_qc_flagged_variants.tsv"), sep = "\t"
)

tail_summary[, label := factor(
  label, levels = rev(quality_metrics$label)
)]
p9 <- ggplot(
  tail_summary,
  aes(fraction_concerning_tail, label, fill = region)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68,
           orientation = "y") +
  geom_vline(xintercept = 0.05, linetype = "dotted",
             linewidth = 0.45, colour = "grey45") +
  facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c("KIR" = "#B03A2E", "Other post_ppi" = "grey70"),
    name = NULL
  ) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = "Variants in the concerning tail of each quality metric",
    subtitle = "Cut-off: 5th or 95th percentile of non-KIR observations; strict boundary, exploratory",
    x = "Percentage of variants", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(image_out, "kir_qc_tail_comparison.png"), p9,
       width = 12, height = 8, dpi = 300)

# Compact gene-level summaries for identifying a few problematic KIR genes.
kir_gene_summary <- qc_sites[IS_KIR == 1L, .(
  n_sites = .N,
  median_MQ = median(MQ, na.rm = TRUE),
  median_QD = median(QD, na.rm = TRUE),
  median_VQSLOD = median(VQSLOD, na.rm = TRUE)
), by = KIR_GENE][order(median_VQSLOD)]
fwrite(kir_gene_summary, file.path(data_in, "kir_qc_by_gene.tsv"), sep = "\t")

# Allele frequency is not a quality metric, but its distribution can confound
# quality comparisons. Describe it separately rather than adjusting silently.
allele_frequency_summary <- qc_sites[, .(
  n_AF_available = sum(is.finite(AF)),
  median_AF = if (any(is.finite(AF))) median(AF[is.finite(AF)]) else NA_real_,
  median_AC = if (any(is.finite(AC))) as.numeric(median(AC[is.finite(AC)])) else NA_real_
), by = region]
fwrite(allele_frequency_summary,
       file.path(data_in, "kir_qc_allele_frequency_summary.tsv"), sep = "\t")

# Produce a readable, local HTML report, with graphics linked from the image folder.
html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}
html_table <- function(x) {
  if (!nrow(x)) return("<p>No observations.</p>")
  cols <- paste0("<th>", html_escape(names(x)), "</th>", collapse = "")
  rows <- vapply(seq_len(nrow(x)), function(i) {
    paste0("<tr>", paste0("<td>", html_escape(unlist(x[i], use.names = FALSE)),
                          "</td>", collapse = ""), "</tr>")
  }, character(1))
  paste0("<table><thead><tr>", cols, "</tr></thead><tbody>",
         paste(rows, collapse = "\n"), "</tbody></table>")
}
fmt <- function(x, digits = 2) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}
fmt_p <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) "NA" else format.pval(v, digits = 2, eps = 1e-4)
  }, character(1))
}
report_stats <- qc_statistics[, .(
  Group = type,
  Metric = label,
  `n KIR` = n_KIR,
  `n other` = n_other,
  `Median KIR` = fmt(median_KIR),
  `Median other` = fmt(median_other),
  `Median difference [95% CI]` = paste0(
    fmt(median_difference), " [", fmt(ci_low), ", ", fmt(ci_high), "]"
  ),
  `Cliff's delta` = fmt(cliffs_delta),
  `BH-adjusted p` = fmt_p(p_FDR)
)]
report_tail <- tail_summary[, .(
  Group = type,
  Metric = as.character(label),
  `Concerning direction` = concerning,
  Region = region,
  `Count in tail` = paste0(n_concerning_tail, "/", n),
  Percentage = paste0(fmt(100 * fraction_concerning_tail, 1), "%"),
  `Reference cut-off` = fmt(cutoff)
)]
report_gene <- kir_gene_summary[, .(
  Gene = KIR_GENE,
  Sites = n_sites,
  `Median MQ` = fmt(median_MQ),
  `Median QD` = fmt(median_QD),
  `Median VQSLOD` = fmt(median_VQSLOD)
)]

site_n_KIR <- qc_sites[IS_KIR == 1L, .N]
site_n_other <- qc_sites[IS_KIR == 0L, .N]
het_n_KIR <- het_sites[IS_KIR == 1L, .N]
het_n_other <- het_sites[IS_KIR == 0L, .N]

# Focused, exploratory tests of the low-quality tails for MQ and VQSLOD.
# The thresholds were defined from the same non-KIR data used as controls.
# Fisher p-values therefore describe a data-dependent exploratory screen, not
# confirmatory tests with externally fixed call-failure thresholds.
wilson_interval <- function(events, total, conf_level = 0.95) {
  if (length(events) != 1L || length(total) != 1L ||
      !is.finite(events) || !is.finite(total) || total < 1L ||
      events < 0 || events > total) {
    stop("Invalid numerator or denominator for Wilson interval.")
  }
  z <- qnorm(1 - (1 - conf_level) / 2)
  p <- events / total
  denom <- 1 + z^2 / total
  centre <- (p + z^2 / (2 * total)) / denom
  half_width <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denom
  c(max(0, centre - half_width), min(1, centre + half_width))
}

# Unpaired Newcombe score interval for a risk difference, based on Wilson limits.
newcombe_difference_interval <- function(k1, n1, k0, n0) {
  p1 <- k1 / n1
  p0 <- k0 / n0
  a <- wilson_interval(k1, n1)
  b <- wilson_interval(k0, n0)
  d <- p1 - p0
  c(
    d - sqrt((p1 - a[1])^2 + (b[2] - p0)^2),
    d + sqrt((a[2] - p1)^2 + (p0 - b[1])^2)
  )
}

key_tail_codes <- c("MQ", "VQSLOD")
key_tail_tests <- rbindlist(lapply(key_tail_codes, function(code_value) {
  x <- tail_summary[type == "Site" & code == code_value]
  k <- x[region == "KIR"]
  o <- x[region == "Other post_ppi"]
  if (nrow(k) != 1L || nrow(o) != 1L)
    stop("Expected exactly two site-level tail groups for ", code_value)
  if (!isTRUE(all.equal(as.numeric(k$cutoff), as.numeric(o$cutoff))))
    stop("Mismatched tail thresholds for ", code_value)
  
  ci_k <- wilson_interval(k$n_concerning_tail, k$n)
  ci_o <- wilson_interval(o$n_concerning_tail, o$n)
  ci_diff <- newcombe_difference_interval(
    k$n_concerning_tail, k$n,
    o$n_concerning_tail, o$n
  )
  counts <- matrix(c(
    k$n_concerning_tail, k$n - k$n_concerning_tail,
    o$n_concerning_tail, o$n - o$n_concerning_tail
  ), nrow = 2L, byrow = TRUE)
  data.table(
    code = code_value,
    label = as.character(k$label),
    cutoff = as.numeric(k$cutoff),
    n_KIR = as.integer(k$n),
    low_KIR = as.integer(k$n_concerning_tail),
    proportion_KIR = as.numeric(k$fraction_concerning_tail),
    KIR_ci_low = ci_k[1L],
    KIR_ci_high = ci_k[2L],
    n_other = as.integer(o$n),
    low_other = as.integer(o$n_concerning_tail),
    proportion_other = as.numeric(o$fraction_concerning_tail),
    other_ci_low = ci_o[1L],
    other_ci_high = ci_o[2L],
    difference = as.numeric(k$fraction_concerning_tail - o$fraction_concerning_tail),
    difference_ci_low = ci_diff[1L],
    difference_ci_high = ci_diff[2L],
    fisher_p = as.numeric(stats::fisher.test(counts, alternative = "two.sided")$p.value)
  )
}))
key_tail_tests[, fisher_p_BH := p.adjust(fisher_p, method = "BH")]
fwrite(key_tail_tests, file.path(data_in, "kir_qc_key_tail_tests.tsv"), sep = "\t")

# Confidence intervals expose the uncertainty hidden by a zero observed count.
key_tail_plot <- rbindlist(list(
  key_tail_tests[, .(
    Metric = label, Region = "KIR", count = low_KIR, n = n_KIR,
    proportion = proportion_KIR,
    ci_low = KIR_ci_low, ci_high = KIR_ci_high
  )],
  key_tail_tests[, .(
    Metric = label, Region = "Other post_ppi", count = low_other, n = n_other,
    proportion = proportion_other,
    ci_low = other_ci_low, ci_high = other_ci_high
  )]
))
key_tail_plot[, Metric := factor(
  Metric, levels = key_tail_tests$label
)]
key_tail_plot[, Region := factor(
  Region, levels = c("KIR", "Other post_ppi")
)]

p10 <- ggplot(key_tail_plot, aes(x = Region, y = proportion, colour = Region)) +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high),
                  linewidth = 0.55, fatten = 2.1) +
  geom_text(aes(label = paste0(count, "/", n)),
            colour = "#30343A", vjust = -0.95, size = 3.3) +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  scale_colour_manual(
    values = c("KIR" = "#B03A2E", "Other post_ppi" = "#777777"),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0.08, 0.22))
  ) +
  labs(
    title = "Low-quality tail frequency among selected variants",
    subtitle = "Points: observed proportion; intervals: 95% Wilson confidence intervals",
    x = NULL, y = "Percentage of variant sites"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(image_out, "kir_qc_key_tail_intervals.png"), p10,
       width = 10.5, height = 5.5, dpi = 300)

# Data-derived, concise report. All test results are recalculated on each run.
report_metric <- function(code_value, type_value = "Site") {
  x <- qc_statistics[type == type_value & code == code_value]
  if (nrow(x) != 1L)
    stop("Expected one statistical result for ", type_value, ": ", code_value)
  x
}
report_tail_metric <- function(code_value) {
  x <- key_tail_tests[code == code_value]
  if (nrow(x) != 1L) stop("Expected one key-tail test for ", code_value)
  x
}
fmt_percent <- function(x, digits = 1L) paste0(fmt(100 * x, digits), "%")
fmt_pp <- function(x, digits = 1L) paste0(fmt(100 * x, digits), " pp")
fmt_p_clause <- function(p, label = "p") {
  if (!is.finite(p)) return(paste0(label, " unavailable"))
  if (p < 1e-4) return(paste0(label, " &lt; 0.0001"))
  paste0(label, " = ", formatC(p, digits = 3, format = "f"))
}
fmt_difference_interval <- function(x) {
  paste0(fmt_pp(x$difference), " [", fmt_pp(x$difference_ci_low),
         ", ", fmt_pp(x$difference_ci_high), "]")
}

mq <- report_metric("MQ")
vqs <- report_metric("VQSLOD")
qd <- report_metric("QD")
fs <- report_metric("FS")
sor <- report_metric("SOR")
qual <- report_metric("QUAL")
dp <- report_metric("DP", "Carrier median")
gq <- report_metric("GQ", "Carrier median")
ab <- report_metric("AB_deviation", "Carrier median")
mq_tail <- report_tail_metric("MQ")
vqs_tail <- report_tail_metric("VQSLOD")

report_key_tails <- key_tail_tests[, .(
  Metric = label,
  `Non-KIR fifth-percentile cut-off` = fmt(cutoff),
  `KIR below cut-off` = paste0(low_KIR, "/", n_KIR, " (", fmt_percent(proportion_KIR), ")"),
  `Other below cut-off` = paste0(low_other, "/", n_other, " (", fmt_percent(proportion_other), ")"),
  `Difference [95% CI]` = vapply(seq_len(.N), function(i) {
    paste0(fmt_pp(difference[i]), " [", fmt_pp(difference_ci_low[i]),
           ", ", fmt_pp(difference_ci_high[i]), "]")
  }, character(1)),
  `Fisher p` = fmt_p(fisher_p),
  `BH-adjusted p` = fmt_p(fisher_p_BH)
)]

# One table for short, useful decisions. A p-value alone is not a quality verdict.
report_overview <- data.table(
  Question = c(
    "Are KIR sites overrepresented in the low-MQ tail?",
    "Are KIR sites overrepresented in the low-VQSLOD tail?",
    "Do the whole-score distributions differ?",
    "Are the other observed QC measures consistently worse?"
  ),
  Observation = c(
    paste0(mq_tail$low_KIR, "/", mq_tail$n_KIR, " versus ",
           mq_tail$low_other, "/", mq_tail$n_other,
           " below MQ ", fmt(mq_tail$cutoff)),
    paste0(vqs_tail$low_KIR, "/", vqs_tail$n_KIR, " versus ",
           vqs_tail$low_other, "/", vqs_tail$n_other,
           " below VQSLOD ", fmt(vqs_tail$cutoff)),
    paste0("MQ ", fmt(mq$median_KIR), " vs ", fmt(mq$median_other),
           "; VQSLOD ", fmt(vqs$median_KIR), " vs ", fmt(vqs$median_other)),
    paste0("QD ", fmt(qd$median_KIR), " vs ", fmt(qd$median_other),
           "; GQ ", fmt(gq$median_KIR), " vs ", fmt(gq$median_other),
           "; |AB - 0.5| ", fmt(ab$median_KIR), " vs ", fmt(ab$median_other))
  ),
  Assessment = c(
    "No observed low-MQ KIR sites; uncertainty remains",
    "Higher observed KIR fraction, with a wide confidence interval",
    "Both MQ and VQSLOD show distributional shifts",
    "No consistent deterioration across these selected calls"
  )
)

summary_rationale <- paste0(
  "Our genome-wide analysis identified enrichment of rare variants in disease-",
  "associated infection-response pathways that included several KIR genes. ",
  "Because KIR genes exhibit extensive sequence homology and copy-number ",
  "variation, we considered whether read-mapping or variant-calling artefacts ",
  "could have contributed to these findings. We therefore performed a post-hoc ",
  "quality assessment of ", site_n_KIR, " KIR-annotated SNVs against ",
  site_n_other, " other retained SNVs from the same three MCL extracts."
)

summary_mapping <- paste0(
  "Median mapping quality was ", fmt(mq$median_KIR), " in KIR versus ",
  fmt(mq$median_other), " in other sites (Wilcoxon ",
  fmt_p_clause(mq$p_FDR, "BH-adjusted p"), "). None of ",
  mq_tail$n_KIR, " KIR sites fell below the other variants' fifth-percentile ",
  "MQ cut-off (", fmt(mq_tail$cutoff), "), versus ",
  mq_tail$low_other, "/", mq_tail$n_other, " non-KIR sites. ",
  "The absolute difference in low-MQ frequency was ",
  fmt_difference_interval(mq_tail), " (Newcombe 95% CI; KIR minus non-KIR). ",
  "The exploratory Fisher test gave ", fmt_p_clause(mq_tail$fisher_p),
  " (two-test BH-adjusted ", fmt_p_clause(mq_tail$fisher_p_BH, "p"), "). ",
  "The observed zero does not imply that mapping errors are absent."
)

summary_recalibration <- paste0(
  "VQSLOD is a model-derived GATK recalibration score combining several quality ",
  "annotations; a negative value may still pass a VQSR tranche ",
  "threshold and does not, by itself, establish an erroneous variant call. ",
  "VQSLOD was lower in KIR (median ", fmt(vqs$median_KIR), " versus ",
  fmt(vqs$median_other), "; Wilcoxon ",
  fmt_p_clause(vqs$p_FDR, "BH-adjusted p"), "). At the other variants' ",
  "fifth-percentile VQSLOD cut-off (", fmt(vqs_tail$cutoff), "), ",
  vqs_tail$low_KIR, "/", vqs_tail$n_KIR, " KIR variants (",
  fmt_percent(vqs_tail$proportion_KIR), ") and ",
  vqs_tail$low_other, "/", vqs_tail$n_other, " other variants (",
  fmt_percent(vqs_tail$proportion_other), ") were in the low-score tail. ",
  "The absolute difference was ", fmt_difference_interval(vqs_tail),
  "; exploratory Fisher ", fmt_p_clause(vqs_tail$fisher_p),
  ". A shift in this recalibration score does not identify individual mapping errors."
)

strand_result <- if (all(is.finite(c(fs$p_FDR, sor$p_FDR))) &&
                     all(c(fs$p_FDR, sor$p_FDR) >= 0.05)) {
  "FS and SOR showed no significant differences after the nine-metric BH correction. "
} else {
  "FS and SOR results are given separately in the full comparison table. "
}

summary_other <- paste0(
  "Quality by depth (median ", fmt(qd$median_KIR), " versus ",
  fmt(qd$median_other), "), per-variant median heterozygous-carrier read ",
  "depth (", fmt(dp$median_KIR), " versus ", fmt(dp$median_other),
  "), median genotype quality (", fmt(gq$median_KIR), " versus ",
  fmt(gq$median_other), ") and median absolute heterozygous allele-balance ",
  "deviation (", fmt(ab$median_KIR), " versus ",
  fmt(ab$median_other), ") did not show a consistent pattern of poor KIR ",
  "quality. ", strand_result, "Variant QUAL was higher in KIR (",
  fmt(qual$median_KIR), " versus ", fmt(qual$median_other),
  "), but QUAL is not independent evidence of correct mapping. ",
  "GQ values of 99 may be capped."
)

pass_note <- if (isTRUE(all(qc_sites$FILTER == "PASS"))) {
  "All retained sites had FILTER=PASS, which cannot discriminate within this selected set. "
} else {
  "Only variants present in the supplied downstream MCL files were evaluated. "
}

summary_conclusion <- paste0(
  "The post-hoc comparison found no consistent evidence of poorer quality ",
  "among the retained KIR calls. Quality by depth (median ",
  fmt(qd$median_KIR), " versus ", fmt(qd$median_other),
  "), strand-bias metrics, genotype quality (median ", fmt(gq$median_KIR),
  " versus ", fmt(gq$median_other), ") and heterozygous allele-balance ",
  "deviation (", fmt(ab$median_KIR), " versus ", fmt(ab$median_other),
  ") were broadly comparable; KIR calls had higher median read depth ",
  "and QUAL scores. Median mapping quality was slightly lower (",
  fmt(mq$median_KIR), " versus ", fmt(mq$median_other),
  "), but ", mq_tail$low_KIR, "/", mq_tail$n_KIR,
  " KIR sites versus ", mq_tail$low_other, "/", mq_tail$n_other,
  " other sites fell below the comparator's fifth-percentile MQ threshold. ",
  "VQSLOD, GATK's model-derived recalibration score, was markedly lower ",
  "(median ", fmt(vqs$median_KIR), " versus ", fmt(vqs$median_other),
  ")", if (is.finite(vqs$p_FDR) && vqs$p_FDR < 0.05)
    " and differed significantly between groups" else "",
  "; however, the proportion in the lowest VQSLOD tail ",
  if (is.finite(vqs_tail$fisher_p) && vqs_tail$fisher_p < 0.05)
    "also differed in this exploratory test" else
      "was not significantly different in this exploratory test",
  " (", fmt_percent(vqs_tail$proportion_KIR), " versus ",
  fmt_percent(vqs_tail$proportion_other), "). ",
  "The available metrics do not indicate a general QC failure that would ",
  "justify excluding KIR variants or changing the pipeline's existing QC ",
  "criteria. The VQSLOD shift warrants attention, and these selected calls ",
  "cannot independently validate mapping accuracy or detection sensitivity ",
  "in the KIR region."
)

# Captions and report sections make the statistical question visible before
# displaying every diagnostic plot. Files remain local relative image links.
html <- paste0(
  '<!doctype html><html lang="en"><head><meta charset="utf-8">',
  '<meta name="viewport" content="width=device-width,initial-scale=1">',
  '<title>KIR variant quality report</title>',
  '<style>',
  ':root{color-scheme:light}',
  'body{font:15px/1.62 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;',
  'color:#202428;background:#fff;max-width:1120px;margin:38px auto;padding:0 24px 70px}',
  'h1{font-size:30px;margin:0 0 8px;letter-spacing:-.02em}',
  'h2{font-size:21px;margin:42px 0 12px;line-height:1.25}',
  'h3{font-size:16px;margin:23px 0 7px;line-height:1.3}',
  'p{max-width:100ch;margin:9px 0 16px}',
  '.note{color:#596169;font-size:13px}',
  '.lead{font-size:16px;max-width:94ch}',
  '.table-wrap{overflow-x:auto;margin:12px 0 22px}',
  'table{border-collapse:collapse;width:100%;font-size:13px}',
  'th,td{padding:9px 12px;border-bottom:1px solid #e1e3e5;text-align:left;vertical-align:top}',
  'th{background:#f5f6f7;font-weight:650;white-space:nowrap}',
  'td{white-space:normal;min-width:75px}',
  '.table-wrap table{margin:0}',
  '.key-figure{margin:23px 0 28px}',
  'figure{margin:24px 0 28px}',
  'figure img{display:block;width:100%;height:auto}',
  'figcaption{font-size:13px;color:#596169;margin-top:8px}',
  'nav{font-size:13px;margin:14px 0 27px}',
  'nav a{color:#314b64;margin-right:16px;text-decoration:none}',
  'nav a:hover{text-decoration:underline}',
  '.separate{border-top:1px solid #e1e3e5;padding-top:4px}',
  '@media print{body{margin:0;padding:0 8mm;max-width:none}',
  'figure,.table-wrap{break-inside:avoid}nav{display:none}}',
  '</style></head><body>',
  '<h1>KIR variant quality assessment</h1>',
  '<p class="note">Exploratory post-filtering QC assessment | MCL_22, MCL_586 and MCL_836 | ',
  site_n_KIR, ' KIR and ', site_n_other, ' other distinct SNVs</p>',
  '<nav><a href="#findings">Findings</a><a href="#tail">Low-score tails</a>',
  '<a href="#distributions">Figures</a><a href="#statistics">Full statistics</a>',
  '<a href="#limitations">Limitations</a></nav>',
  
  '<section id="findings">',
  '<h2>Research question and conclusion</h2>',
  '<p>', html_escape(summary_rationale), '</p>',
  '<p class="lead"><strong>Conclusion.</strong> ', html_escape(summary_conclusion), '</p>',
  '<h2>Principal results</h2>',
  '<div class="table-wrap">', html_table(report_overview), '</div>',
  '<h3>Mapping quality</h3>',
  '<p>', summary_mapping, '</p>',
  '<h3>Recalibrated variant quality</h3>',
  '<p>', summary_recalibration, '</p>',
  '<h3>Other measures</h3>',
  '<p>', summary_other, '</p>',
  '</section>',
  
  '<section id="tail" class="separate">',
  '<h2>Low-score tails: frequency and uncertainty</h2>',
  '<p>The original concern concerns an excess of low-quality KIR calls. ',
  'For this question, the frequency of observations in the concerning tail ',
  'is more interpretable than comparing medians alone. These are ',
  '<strong>relative</strong> cut-offs derived from the non-KIR sample, ',
  'not established mapping-failure thresholds.</p>',
  '<figure class="key-figure"><img src="kir_qc_key_tail_intervals.png" ',
  'alt="Low-MQ and low-VQSLOD proportions with Wilson confidence intervals">',
  '<figcaption>Site-level rates and two-sided 95% Wilson intervals. ',
  'Each SNV contributes once; a zero observed count still has an upper ',
  'confidence limit.</figcaption></figure>',
  '<div class="table-wrap">', html_table(report_key_tails), '</div>',
  '<p class="note">Risk differences are KIR minus other variants in ',
  'percentage points, with 95% Newcombe score confidence intervals. ',
  'Fisher tests are two-sided; BH correction covers these two exploratory ',
  'tail comparisons only. The threshold was estimated using the same ',
  'non-KIR sample being tested, so the nominal p-values should not be ',
  'interpreted as confirmatory error rates.</p>',
  '</section>',
  
  '<section id="distributions" class="separate">',
  '<h2>Quality distributions and individual variants</h2>',
  '<p>The complete distributions distinguish shifts affecting many variants ',
  'from a small number of particularly low-scoring observations. ',
  'KIR genes are coloured individually; non-KIR variants appear in grey.</p>',
  '<figure><img src="kir_site_quality_histograms.png" ',
  'alt="Site-quality distributions with KIR variants marked">',
  '<figcaption>Non-KIR variant counts in grey; individual KIR sites marked ',
  'by gene along the horizontal axis.</figcaption></figure>',
  '<figure><img src="kir_quality_space.png" ',
  'alt="Quality by depth versus VQSLOD scatter plot">',
  '<figcaption>Individual SNVs in QD versus VQSLOD space.</figcaption></figure>',
  '<figure><img src="kir_chr19_quality.png" ',
  'alt="Chromosome 19 quality at individual SNVs">',
  '<figcaption>Only the MCL extracts are used. All observed chromosome 19 ',
  'variants are points; no reference VCF, smoothed line or window median is included.</figcaption></figure>',
  '<figure><img src="kir_chr19_quality_zoom.png" ',
  'alt="Local chromosome 19 KIR quality">',
  '<figcaption>Chromosome 19 view enlarged around the KIR sites.</figcaption></figure>',
  '<figure><img src="kir_all_site_quality_ranked.png" ',
  'alt="Quality scores ranked among retained variants">',
  '<figcaption>Percentiles are computed separately for each metric.</figcaption></figure>',
  '<h3>Heterozygous carrier quality</h3>',
  '<figure><img src="kir_carrier_quality_histograms.png" ',
  'alt="Per-variant median heterozygous carrier-quality histograms">',
  '<figcaption>Each distinct SNV has one median of its heterozygous ',
  'carrier calls, avoiding extra weight for variants seen in more samples.</figcaption></figure>',
  '<figure><img src="kir_carrier_quality_ranked.png" ',
  'alt="Per-variant median carrier quality ranks"></figure>',
  '</section>',
  
  '<section id="statistics" class="separate">',
  '<h2>Detailed statistical results</h2>',
  '<h3>Descriptive summaries</h3>',
  '<p>Distinct variant sites:</p>',
  '<div class="table-wrap">', html_table(site_summary), '</div>',
  '<p>Heterozygous genotype observations (not independent for inference):</p>',
  '<div class="table-wrap">', html_table(carrier_summary), '</div>',
  '<h3>Complete score-distribution comparisons</h3>',
  '<p>Site-level analyses use distinct SNVs. Carrier-level analyses use ',
  'one median heterozygous genotype per distinct SNV. ',
  'Two-sided Wilcoxon rank-sum tests were adjusted together across ',
  nrow(qc_statistics), ' metrics using Benjamini-Hochberg. ',
  'The confidence interval is for the difference in medians ',
  '(KIR minus other) from 1,500 variant-level bootstrap resamples. ',
  'Cliff&#39;s delta measures relative rank separation; it is not a ',
  'variant-failure probability.</p>',
  '<div class="table-wrap">', html_table(report_stats), '</div>',
  '<h3>Additional relative-tail screens</h3>',
  '<p>For each metric, the concerning direction was defined as below ',
  'the fifth percentile (MQ, QD, VQSLOD, QUAL, DP and GQ) or above the ',
  '95th percentile (FS, SOR and allele-balance deviation) of the ',
  'non-KIR values. Ties and strict inequalities can change the baseline ',
  'fraction from exactly 5%.</p>',
  '<figure><img src="kir_qc_tail_comparison.png" ',
  'alt="Exploratory concerning-tail rates for all QC metrics"></figure>',
  '<div class="table-wrap">', html_table(report_tail), '</div>',
  '<h3>Breakdown by KIR gene</h3>',
  '<p>Sites annotated to several KIR genes form a separate group, ',
  'rather than being assigned arbitrarily.</p>',
  '<div class="table-wrap">', html_table(report_gene), '</div>',
  '<h3>Allele-frequency composition</h3>',
  '<p>AF and AC may confound comparison of the quality distributions ',
  'but are not themselves measurements of read quality.</p>',
  '<div class="table-wrap">', html_table(allele_frequency_summary), '</div>',
  '</section>',
  
  '<section id="limitations" class="separate">',
  '<h2>Interpretation and limitations</h2>',
  '<p><strong>Selected variants.</strong> Only biallelic SNVs present in ',
  'the three downstream MCL extracts are assessed. The complete ',
  'chromosome 19 call set, non-PASS variants, missing calls and false ',
  'negatives are not represented. ', pass_note, '</p>',
  '<p><strong>Genomic structure.</strong> KIR paralogy, copy-number ',
  'variation, reference representation and read placement can produce ',
  'confident but incorrect alignments. The available VCF metrics cannot ',
  'independently establish mapping correctness. Nearby variants may ',
  'not be statistically independent, so site-resampling intervals and ',
  'rank tests may underestimate uncertainty.</p>',
  '<p><strong>Recalibration.</strong> VQSLOD is dependent on the original ',
  'VQSR training model and annotations; neither its negative sign nor ',
  'the empirical fifth-percentile cut-off automatically identifies a ',
  'false positive. The model, culprit annotations and BAM/CRAM evidence ',
  'were not analysed. QUAL, GQ and DP are not independent confirmations ',
  'of correct gene assignment.</p>',
  '<p><strong>Inference.</strong> Non-significant tests do not prove ',
  'equivalence. These tail thresholds and analyses are exploratory, ',
  'and the variants are not a random sample of all callable exome sites. ',
  'Genomic region and allele-frequency composition may differ between ',
  'groups. A confirmatory QC claim would require independent validation ',
  'or a pre-specified, externally justified QC threshold.</p>',
  '<p class="note">Supporting TSV outputs: ', html_escape(data_in),
  '. Figures and report: ', html_escape(image_out), '.</p>',
  '</section></body></html>'
)

writeLines(html, file.path(image_out, "kir_qc_report.html"), useBytes = TRUE)
cat("\nExploratory MQ and VQSLOD low-tail comparisons:\n")
print(key_tail_tests[, .(
  code, low_KIR, n_KIR, low_other, n_other,
  difference, difference_ci_low, difference_ci_high, fisher_p, fisher_p_BH
)])
cat("\nQC statistical comparison (post_ppi):\n")
print(qc_statistics[, .(
  type, code, n_KIR, n_other, median_KIR, median_other,
  median_difference, ci_low, ci_high, cliffs_delta, p_value, p_FDR
)])
cat("\nWrote report: ", file.path(image_out, "kir_qc_report.html"), "\n", sep = "")
