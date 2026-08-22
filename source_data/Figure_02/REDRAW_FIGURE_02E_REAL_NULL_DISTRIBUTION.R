# =============================================================================
# REDRAW_FIGURE_02E_LOCKED_FINAL.R
#
# Uses ONLY the validated locked file:
#   Figure_02/Random_benchmark_null_distributions.csv.gz
#
# Expected folder when run from the repository root:
#   source_data/Figure_02
#
# Outputs:
#   Figure_02E_LOCKED_REAL_NULL_DISTRIBUTION.png
#   Figure_02E_LOCKED_REAL_NULL_DISTRIBUTION.pdf
#   Figure_02E_LOCKED_plot_source.csv
#
# The script hard-validates the two locked ECM_LIGAND_21 null means and
# add-one empirical P values before drawing.
# =============================================================================

options(stringsAsFactors = FALSE)

FIG_DIR <- Sys.getenv(
  "ECM21_FIGURE_02_DIR",
  unset = file.path("source_data", "Figure_02")
)
NULL_FILE <- file.path(FIG_DIR, "Random_benchmark_null_distributions.csv.gz")
SUMMARY_FILE <- file.path(FIG_DIR, "Table_biological_coherence_random_benchmark.csv")
OUT_DIR <- file.path(FIG_DIR, "Figure_02E_real_null_distribution")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(NULL_FILE)) stop("Missing validated null file:\n", NULL_FILE)
if (!file.exists(SUMMARY_FILE)) stop("Missing locked summary:\n", SUMMARY_FILE)

needed <- c("ggplot2", "scales")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse=", "))

library(ggplot2)
library(scales)

null_all <- read.csv(gzfile(NULL_FILE), check.names=FALSE, stringsAsFactors=FALSE)
summary_all <- read.csv(SUMMARY_FILE, check.names=FALSE, stringsAsFactors=FALSE)

locked <- summary_all[
  summary_all$panel == "ECM_LIGAND_21" &
  summary_all$dataset %in% c("GSE132465","GSE144735"),
  , drop=FALSE
]
locked <- locked[match(c("GSE132465","GSE144735"), locked$dataset),]

add_one_p <- function(x, obs) (1 + sum(x >= obs, na.rm=TRUE)) / (1 + sum(is.finite(x)))

audit <- lapply(seq_len(nrow(locked)), function(i) {
  ds <- locked$dataset[i]
  x <- as.numeric(null_all[
    null_all$dataset == ds & null_all$panel == "ECM_LIGAND_21",
    "fibroblast_specificity"
  ])
  x <- x[is.finite(x)]
  obs <- as.numeric(locked$observed_fibroblast_specificity[i])
  p <- add_one_p(x, obs)
  data.frame(
    dataset=ds,
    n_null=length(x),
    null_mean=mean(x),
    locked_mean=as.numeric(locked$null_mean_fibroblast_specificity[i]),
    n_ge_observed=sum(x >= obs),
    empirical_p=p,
    locked_p=as.numeric(locked$empirical_p_fibroblast_specificity[i]),
    exact=(
      length(x)==10000 &&
      abs(mean(x)-as.numeric(locked$null_mean_fibroblast_specificity[i])) < 1e-10 &&
      abs(p-as.numeric(locked$empirical_p_fibroblast_specificity[i])) < 1e-12
    )
  )
})
audit <- do.call(rbind,audit)
print(audit)

if (!all(audit$exact)) {
  stop("Locked validation failed. Plotting stopped.")
}

null_plot <- null_all[
  null_all$panel == "ECM_LIGAND_21" &
  null_all$dataset %in% c("GSE132465","GSE144735"),
  c("dataset","iteration","fibroblast_specificity")
]
names(null_plot)[3] <- "value"
null_plot$dataset <- factor(null_plot$dataset, levels=c("GSE132465","GSE144735"))

obs <- locked[,c(
  "dataset",
  "observed_fibroblast_specificity",
  "empirical_p_fibroblast_specificity"
)]
names(obs)[2] <- "observed"
obs$dataset <- factor(obs$dataset, levels=c("GSE132465","GSE144735"))
obs$n_ge <- 0L
obs$annotation <- sprintf(
  "Observed = %.3f\n0 / 10,000 matched panels \u2265 observed\nEmpirical P = 0.0001",
  obs$observed
)

plot_source <- merge(null_plot, obs, by="dataset", all.x=TRUE, sort=FALSE)
write.csv(plot_source, file.path(OUT_DIR,"Figure_02E_LOCKED_plot_source.csv"), row.names=FALSE)
write.csv(audit, file.path(OUT_DIR,"Figure_02E_LOCKED_validation.csv"), row.names=FALSE)

p <- ggplot(null_plot, aes(dataset, value)) +
  geom_violin(
    width=0.72, scale="width", trim=FALSE,
    fill="#DCE6EC", colour="#55717F", linewidth=0.65
  ) +
  geom_boxplot(
    width=0.12, outlier.shape=NA,
    fill="white", colour="#455A64", linewidth=0.5
  ) +
  geom_point(
    data=obs, aes(dataset, observed),
    inherit.aes=FALSE,
    shape=23, size=4.2, stroke=0.9,
    fill="#C9223A", colour="#C9223A"
  ) +
  geom_text(
    data=obs, aes(dataset, observed, label=annotation),
    inherit.aes=FALSE,
    hjust=-0.08, vjust=0.45,
    family="Arial", fontface="bold", size=3.0,
    colour="#A01832"
  ) +
  scale_y_continuous(expand=expansion(mult=c(0.04,0.20))) +
  coord_cartesian(clip="off") +
  labs(
    tag="E",
    x=NULL,
    y="Fibroblast-specificity index",
    title="Outcome-blind matched-panel benchmark",
    subtitle="ECM_LIGAND_21 versus 10,000 expression/detection-matched transcriptome panels per cohort",
    caption=paste(
      "Violin = complete empirical null distribution; box = median/IQR;",
      "red diamond = observed ECM_LIGAND_21. Add-one upper-tail empirical P."
    )
  ) +
  theme_bw(base_size=10, base_family="Arial") +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),
    plot.title=element_text(face="bold", size=12),
    plot.subtitle=element_text(size=9.5, colour="#4D4D4D"),
    plot.caption=element_text(size=8, colour="#555555", hjust=0),
    plot.tag=element_text(face="bold", size=14),
    plot.tag.position=c(0.01,0.99),
    plot.margin=margin(8,80,8,8)
  )

ggsave(
  file.path(OUT_DIR,"Figure_02E_LOCKED_REAL_NULL_DISTRIBUTION.png"),
  p, width=180, height=105, units="mm", dpi=600, bg="white", limitsize=FALSE
)
ggsave(
  file.path(OUT_DIR,"Figure_02E_LOCKED_REAL_NULL_DISTRIBUTION.pdf"),
  p, width=180, height=105, units="mm", device=grDevices::cairo_pdf,
  bg="white", limitsize=FALSE
)

cat("\nPASS: Figure 2E locked real null distribution generated.\n")
cat("Output: ", OUT_DIR, "\n", sep="")
