if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "rtracklayer",
               "eulerr",
               "ggrastr",
               "DESeq2",
               "ggpubr",
               "wigglescout",
               "bedscout"
)

# disable scientific notation
options(scipen=999)

# export 
result_folder = "../results/wigglescout/"

bws <- list.files("../data/CutNTag_ChIP-Seq/bw/",
                  full.names = TRUE)

# diff signal function
bw_granges_diff_analysis <- function(granges_c1,
                                     granges_c2,
                                     label_c1,
                                     label_c2,
                                     estimate_size_factors = FALSE,
                                     as_granges = FALSE) {
  # Bind first, get numbers after
  names_values <- NULL
  fields <- names(mcols(granges_c1))
  
  if ("name" %in% fields) {
    names_values <- mcols(granges_c1)[["name"]]
    granges_c1 <- granges_c1[, fields[fields != "name"]]
  }
  
  fields <- names(mcols(granges_c2))
  if ("name" %in% fields) {
    granges_c2 <- granges_c2[, fields[fields != "name"]]
  }
  
  cts_df <- cbind(data.frame(granges_c1), mcols(granges_c2))
  
  if (!is.null(names_values)) {
    rownames(cts_df) <- names_values
  }
  
  # Needs to drop non-complete cases and match rows
  complete <- complete.cases(cts_df)
  cts_df <- cts_df[complete, ]
  
  values_df <- cts_df[, 6:ncol(cts_df)] %>% dplyr::select(where(is.numeric))
  cts <- get_nreads_columns(values_df)
  
  condition_labels <- c(rep(label_c1, length(mcols(granges_c1))), rep(label_c2, length(mcols(granges_c2))))
  
  
  coldata <- data.frame(colnames(cts), condition = as.factor(condition_labels))
  print(coldata)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~ condition,
    rowRanges = granges_c1[complete, ]
  )
  
  
  if (estimate_size_factors == TRUE) {
    dds <- DESeq2::estimateSizeFactors(dds)
  }
  else {
    # Since files are scaled, we do not want to estimate size factors
    sizeFactors(dds) <- c(rep(1, ncol(cts)))
  }
  
  dds <- DESeq2::estimateDispersions(dds)
  dds <- DESeq2::nbinomWaldTest(dds)
  
  if (as_granges) {
    result <- DESeq2::results(dds, format = "GRanges", alpha = 0.01)
    if (!is.null(names_values)) {
      result$name <- names_values[complete]
    }
    
  }
  else {
    # result <- results(dds, format="DataFrame")
    result <- dds
  }
  
  result
}

get_nreads_columns <- function(df, length_factor = 100) {
  # Convert mean coverages to round integer read numbers
  cts <- as.matrix(df)
  cts <- as.matrix(cts[complete.cases(cts), ])
  cts <- round(cts * length_factor)
  cts
}

compute_logfc = function(
    control_bws,
    control_label,
    trt_bws,
    trt_label,
    peak_subset,
    peak_label) {
  
  # control
  control <- bw_loci(control_bws, peak_subset)
  control$name = paste("control", as.character(seq(1:(length(
    control
  )))), sep = "_")
  
  # sample
  trt <- bw_loci(trt_bws, peak_subset)
  trt$name = paste("treatment", as.character(seq(1:(length(
    trt
  )))), sep = "_")
  
  diff <- bw_granges_diff_analysis(granges_c1 = trt, granges_c2 = control, 
                                   "trt", "control")
  
  lfc = DESeq2::lfcShrink(diff, coef = "condition_trt_vs_control", type =
                            "apeglm")
  fcs = as_tibble(lfc) # keep DESeq2 output table
  return(fcs)
  
}

compute_aggr_logfc = function(
    control_bws,
    control_label,
    trt_bws,
    trt_label,
    peak_subset,
    peak_label) {
  
  # control
  control <- bw_loci(control_bws, peak_subset)
  control$name = paste("control", as.character(seq(1:(length(
    control
  )))), sep = "_")
  
  # sample
  trt <- bw_loci(trt_bws, peak_subset)
  trt$name = paste("treatment", as.character(seq(1:(length(
    trt
  )))), sep = "_")
  
  diff <- bw_granges_diff_analysis(granges_c1 = trt, granges_c2 = control, 
                                               "trt", "control")
  
  lfc = DESeq2::lfcShrink(diff, coef = "condition_trt_vs_control", type =
                                         "apeglm")
  fcs = as_tibble(lfc) %>% filter(pvalue < 0.05) # keep only significant differences
  fcs = fcs$log2FoldChange
  mean_fc = mean(fcs[!is.na(fcs)])
  sd_fc = sd(fcs[!is.na(fcs)])
  
  output = tibble(mean_log2FC = mean_fc, sd_log2FC = sd_fc, DESeq2_contrast = 
                    paste0(trt_label, "_to_", control_label), subset = peak_label)
  
  return(output)
  
}

# Wulfridge CTCF ChIP-Seq mocks
mocks = c("../data/CutNTag_ChIP-Seq/bw/SRR23958386_GSM7116278_E14_Mock_CTCF_Rep2_Mus_musculus_ChIP-Seq.CPMnorm.bw",
          "../data/CutNTag_ChIP-Seq/bw/SRR23958387_GSM7116277_E14_Mock_CTCF_Rep1_Mus_musculus_ChIP-Seq.CPMnorm.bw")

# PDS
# Wulfridge CTCF ChIP-Seq
pds_bigwigs = c("../data/CutNTag_ChIP-Seq/bw/SRR23958384_GSM7116280_E14_PDS_CTCF_Rep2_Mus_musculus_ChIP-Seq.CPMnorm.bw",
                "../data/CutNTag_ChIP-Seq/bw/SRR23958385_GSM7116279_E14_PDS_CTCF_Rep1_Mus_musculus_ChIP-Seq.CPMnorm.bw")

### Compute mean aggregated logFC tables
pds_proms_CTCFG4_aggr_log2fc = compute_aggr_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = pds_bigwigs,
  trt_label = "PDS",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_with_promoters.sorted-CTCF_G4.bed",
  peak_label = "prom_CTCF_G4"
)
pds_proms_CTCFonly_aggr_log2fc = compute_aggr_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = pds_bigwigs,
  trt_label = "PDS",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_with_promoters.sorted-CTCFonly.bed",
  peak_label = "prom_CTCFonly"
)
pds_wo_proms_CTCFG4_aggr_log2fc = compute_aggr_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = pds_bigwigs,
  trt_label = "PDS",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_without_promoters.sorted-CTCF_G4.bed",
  peak_label = "without_prom_CTCF_G4"
)
pds_wo_proms_CTCFonly_aggr_log2fc = compute_aggr_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = pds_bigwigs,
  trt_label = "PDS",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_without_promoters.sorted-CTCFonly.bed",
  peak_label = "without_prom_CTCFonly"
)

pds_log2fc_aggr = rbind(pds_proms_CTCFG4_aggr_log2fc, pds_proms_CTCFonly_aggr_log2fc, pds_wo_proms_CTCFG4_aggr_log2fc,
                   pds_wo_proms_CTCFonly_aggr_log2fc)
# PhenDC3
# Wulfridge CTCF ChIP-Seq
phendc3_bigwigs = c("../data/CutNTag_ChIP-Seq/bw/SRR23958382_GSM7116282_E14_PhenDC3_CTCF_Rep2_Mus_musculus_ChIP-Seq.CPMnorm.bw",
                    "../data/CutNTag_ChIP-Seq/bw/SRR23958383_GSM7116281_E14_PhenDC3_CTCF_Rep1_Mus_musculus_ChIP-Seq.CPMnorm.bw")

phendc3_proms_CTCFG4_aggr_log2fc = compute_aggr_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = phendc3_bigwigs,
  trt_label = "PhenDC3",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_with_promoters.sorted-CTCF_G4.bed",
  peak_label = "prom_CTCF_G4"
)
phendc3_proms_CTCFonly_aggr_log2fc = compute_aggr_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = phendc3_bigwigs,
  trt_label = "PhenDC3",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_with_promoters.sorted-CTCFonly.bed",
  peak_label = "prom_CTCFonly"
)
phendc3_wo_proms_CTCFG4_aggr_log2fc = compute_aggr_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = phendc3_bigwigs,
  trt_label = "PhenDC3",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_without_promoters.sorted-CTCF_G4.bed",
  peak_label = "without_prom_CTCF_G4"
)
phendc3_wo_proms_CTCFonly_aggr_log2fc = compute_aggr_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = phendc3_bigwigs,
  trt_label = "PhenDC3",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_without_promoters.sorted-CTCFonly.bed",
  peak_label = "without_prom_CTCFonly"
)

phendc3_log2fc_aggr = rbind(phendc3_proms_CTCFG4_aggr_log2fc, phendc3_proms_CTCFonly_aggr_log2fc, phendc3_wo_proms_CTCFG4_aggr_log2fc,
                       phendc3_wo_proms_CTCFonly_aggr_log2fc)

wulfridge_log2fc = rbind(pds_log2fc_aggr, phendc3_log2fc_aggr) %>% 
  mutate(source = "Wulfridge_et_al-CTCF_ChIPSeq")

### Compute logFC tables
# PDS
pds_proms_CTCFG4_log2fc = compute_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = pds_bigwigs,
  trt_label = "PDS",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_with_promoters.sorted-CTCF_G4.bed",
  peak_label = "prom_CTCF_G4"
)
pds_proms_CTCFonly_log2fc = compute_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = pds_bigwigs,
  trt_label = "PDS",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_with_promoters.sorted-CTCFonly.bed",
  peak_label = "prom_CTCFonly"
)
pds_wo_proms_CTCFG4_log2fc = compute_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = pds_bigwigs,
  trt_label = "PDS",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_without_promoters.sorted-CTCF_G4.bed",
  peak_label = "without_prom_CTCF_G4"
)
pds_wo_proms_CTCFonly_log2fc = compute_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = pds_bigwigs,
  trt_label = "PDS",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_without_promoters.sorted-CTCFonly.bed",
  peak_label = "without_prom_CTCFonly"
)

# PhenDC3
phendc3_proms_CTCFG4_log2fc = compute_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = phendc3_bigwigs,
  trt_label = "PhenDC3",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_with_promoters.sorted-CTCF_G4.bed",
  peak_label = "prom_CTCF_G4"
)
phendc3_proms_CTCFonly_log2fc = compute_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = phendc3_bigwigs,
  trt_label = "PhenDC3",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_with_promoters.sorted-CTCFonly.bed",
  peak_label = "prom_CTCFonly"
)
phendc3_wo_proms_CTCFG4_log2fc = compute_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = phendc3_bigwigs,
  trt_label = "PhenDC3",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_without_promoters.sorted-CTCF_G4.bed",
  peak_label = "without_prom_CTCF_G4"
)
phendc3_wo_proms_CTCFonly_log2fc = compute_logfc(
  control_bws = mocks,
  control_label = "mock",
  trt_bws = phendc3_bigwigs,
  trt_label = "PhenDC3",
  peak_subset =
    "../data/CutNTag_ChIP-Seq/bed/Wulfridge_CTCF_in_6_categories_without_promoters.sorted-CTCFonly.bed",
  peak_label = "without_prom_CTCFonly"
)
