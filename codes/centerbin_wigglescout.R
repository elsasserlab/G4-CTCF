suppressPackageStartupMessages({
  library("data.table")
  library("tidyverse")
  library("ggplot2")
  library("ggpubr")
  library("glue")
  library("wigglescout")
  library("GenomicRanges")
})

mypal = c(rep("cornflowerblue", 1), rep("orange", 1), rep("red2", 1))

# bigwigs + peaks
bws = list.files(path = "../data/CutNTag_ChIP-Seq/bw/", pattern = "SRR*", full.names = TRUE)
peaks = list.files(path = "../data/CutNTag_ChIP-Seq/bed/", pattern = "_sorted", full.names = TRUE)

# helper functions
# define centerbin and neighbouring regions (bin argument)
around_centerbin = function(peak,
                            column_chr = "V1",
                            column1 = "V2",
                            column2 = "V3",
                            bin) {
  
  peak = fread(peak)
  
  peak$length = peak[[column2]] - peak[[column1]]
  print(glue("Average peak length: {mean(peak$length)}"))
  peak$chr = peak[[column_chr]]
  
  peak = peak %>%
    mutate(around_peakcenter_minus = (peak[[column1]] + length / 2) - bin) %>%
    mutate(around_peakcenter_plus = (peak[[column1]] + length / 2) + bin)
  
  peak$name = glue("plusminus_{bin}")
  peak = GRanges(
    seqnames = peak[[column_chr]],
    ranges = IRanges(
      start = peak$around_peakcenter_minus,
      end = peak$around_peakcenter_plus,
      names = peak$name
    )
  )
  print(glue("{unique(peak@ranges@width)} bp around centerbin"))
  return(peak)
}

# bw_loci function for signal aggregation
aggr_signal = function(bw,
                       peak,
                       peak_name,
                       bw_name,
                       bin) {
  aggr = bw_loci(bwfiles = bw, loci = around_centerbin(peak = peak, bin = bin))
  aggr = as_tibble(aggr)
  # colum 6 contains the mean aggregated signal
  aggr = tibble(aggr_signal = aggr[[6]], peak_set = peak_name, sample_name = bw_name)
  return(aggr)
}

### for looping over the promoter peaks and bigwigs
prom_peaks = peaks[grep("with_promoters", peaks)]
peak_names = c("CTCF-only", "G4-only", "G4_CTCF")
bw_rep1 = bws[grep("_Rep1_", bws)]
bw_rep1_name = c("PhenDC3_rep1", "PDS_rep1", "mock_rep1")
bw_rep2 = bws[grep("_Rep2_", bws)]
bw_rep2_name = c("PhenDC3_rep2", "PDS_rep2", "mock_rep2")

### collect signals around peak centers using the promoter specific peak sets
# 20 bp window
centerbin_plusminus_20_rep1 = list()
for(i in seq(bw_rep1)) {
  for(j in seq(prom_peaks)) {
    df = aggr_signal(bw = bw_rep1[i], peak = prom_peaks[j], bw_name = bw_rep1_name[i],
                     peak_name = peak_names[j], bin = 10)
    centerbin_plusminus_20_rep1[[paste0(peak_names[j], "_", bw_rep1_name[i])]] = df
  }
}
centerbin_plusminus_20_rep1 = bind_rows(centerbin_plusminus_20_rep1)
centerbin_plusminus_20_rep1 = centerbin_plusminus_20_rep1 %>% 
  mutate(window_width = "20") %>% 
  mutate(rep = "rep1") 

centerbin_plusminus_20_rep2 = list()
for(i in seq(bw_rep2)) {
  for(j in seq(prom_peaks)) {
    df = aggr_signal(bw = bw_rep2[i], peak = prom_peaks[j], bw_name = bw_rep2_name[i],
                     peak_name = peak_names[j], bin = 10)
    centerbin_plusminus_20_rep2[[paste0(peak_names[j], "_", bw_rep2_name[i])]] = df
  }
}
centerbin_plusminus_20_rep2 = bind_rows(centerbin_plusminus_20_rep2)
centerbin_plusminus_20_rep2 = centerbin_plusminus_20_rep2 %>% 
  mutate(window_width = "20") %>% 
  mutate(rep = "rep2") 

centerbin_plusminus_20 = rbind(centerbin_plusminus_20_rep1, centerbin_plusminus_20_rep2)

# 40 bp window
centerbin_plusminus_40_rep1 = list()
for(i in seq(bw_rep1)) {
  for(j in seq(prom_peaks)) {
    df = aggr_signal(bw = bw_rep1[i], peak = prom_peaks[j], bw_name = bw_rep1_name[i],
                     peak_name = peak_names[j], bin = 20)
    centerbin_plusminus_40_rep1[[paste0(peak_names[j], "_", bw_rep1_name[i])]] = df
  }
}
centerbin_plusminus_40_rep1 = bind_rows(centerbin_plusminus_40_rep1)
centerbin_plusminus_40_rep1 = centerbin_plusminus_40_rep1 %>% 
  mutate(window_width = "40") %>% 
  mutate(rep = "rep1") 

centerbin_plusminus_40_rep2 = list()
for(i in seq(bw_rep2)) {
  for(j in seq(prom_peaks)) {
    df = aggr_signal(bw = bw_rep2[i], peak = prom_peaks[j], bw_name = bw_rep2_name[i],
                     peak_name = peak_names[j], bin = 20)
    centerbin_plusminus_40_rep2[[paste0(peak_names[j], "_", bw_rep2_name[i])]] = df
  }
}
centerbin_plusminus_40_rep2 = bind_rows(centerbin_plusminus_40_rep2)
centerbin_plusminus_40_rep2 = centerbin_plusminus_40_rep2 %>% 
  mutate(window_width = "40") %>% 
  mutate(rep = "rep2") 

centerbin_plusminus_40 = rbind(centerbin_plusminus_40_rep1, centerbin_plusminus_40_rep2)

# 100 bp window
centerbin_plusminus_100_rep1 = list()
for(i in seq(bw_rep1)) {
  for(j in seq(prom_peaks)) {
    df = aggr_signal(bw = bw_rep1[i], peak = prom_peaks[j], bw_name = bw_rep1_name[i],
                     peak_name = peak_names[j], bin = 50)
    centerbin_plusminus_100_rep1[[paste0(peak_names[j], "_", bw_rep1_name[i])]] = df
  }
}
centerbin_plusminus_100_rep1 = bind_rows(centerbin_plusminus_100_rep1)
centerbin_plusminus_100_rep1 = centerbin_plusminus_100_rep1 %>% 
  mutate(window_width = "100") %>% 
  mutate(rep = "rep1") 

centerbin_plusminus_100_rep2 = list()
for(i in seq(bw_rep2)) {
  for(j in seq(prom_peaks)) {
    df = aggr_signal(bw = bw_rep2[i], peak = prom_peaks[j], bw_name = bw_rep2_name[i],
                     peak_name = peak_names[j], bin = 50)
    centerbin_plusminus_100_rep2[[paste0(peak_names[j], "_", bw_rep2_name[i])]] = df
  }
}
centerbin_plusminus_100_rep2 = bind_rows(centerbin_plusminus_100_rep2)
centerbin_plusminus_100_rep2 = centerbin_plusminus_100_rep2 %>% 
  mutate(window_width = "100") %>% 
  mutate(rep = "rep2") 

centerbin_plusminus_100 = rbind(centerbin_plusminus_100_rep1, centerbin_plusminus_100_rep2)

# combine all signal aggregation
centerbin_plusminus = rbind(centerbin_plusminus_20,centerbin_plusminus_40, centerbin_plusminus_100)
centerbin_plusminus = centerbin_plusminus %>% 
  mutate(peak_type = "with_promoters") %>% 
  separate(sample_name, into = c("sample", "replicate"), 
           sep = "_", remove = FALSE)
rm(list=ls(pattern="^centerbin_plusminus_"))

write_tsv(centerbin_plusminus, "../results/cutntag/centerbins-ws_mean-aggr-promoter_spec_peaks.tsv")

# vis
p1 <- ggviolin(
  centerbin_plusminus,
  x = "window_width",
  y = "aggr_signal",
  fill = "sample",
  palette = mypal,
  add = "median_iqr"
) +
  ylim(0, 10) +
  labs(
    title = "",
    x = "bin size around centerbin",
    y = "mean aggregated signal",
    fill = " "
  ) +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) 
print(p1)

ggsave(
  plot = p1,
  filename = "../results/cutntag/centerbins-ws_mean-aggr-promoter_spec_peaks.pdf",
  device = "pdf",
  width = 7,
  height = 7
)

mypal = c("cornflowerblue", "lightblue", "orange", "#fff7bc", "red2", "#fc9272")
p2 <- ggviolin(
  centerbin_plusminus,
  x = "window_width",
  y = "aggr_signal",
  fill = "sample_name",
  palette = mypal,
  add = "median_iqr"
) +
  ylim(0, 10) +
  labs(
    title = "",
    x = "bin size around centerbin",
    y = "mean aggregated signal",
    fill = " "
  ) +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) 
print(p2)

ggsave(
  plot = p2,
  filename = "../results/cutntag/centerbins-ws_mean-aggr-promoter_spec_peaks-repl_level.pdf",
  device = "pdf",
  width = 7,
  height = 7
)
