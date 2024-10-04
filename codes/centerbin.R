suppressPackageStartupMessages({
  library("data.table")
  library("tidyverse")
  library("ggplot2")
  library("ggpubr")
  library("glue")
})

mypal = c("cornflowerblue", "orange", "red2", "darkgreen", "#505050")

# function returning centerbins based on the deeptools matrix (gz format)
# ideal sample number: 3
# ideal order: MOCK, PDS, PHENDC3
extract_centerbin = function(deeptools_matrix, plot = TRUE) {
  if (ncol(deeptools_matrix) != 1206) {
    stop(
      "The number of columns are not compatible with the Mock/PDS/PhenDC configuration.
          Optimal column number is 1206"
    )
  }
  # if (sum(str_detect(deeptools_matrix$V4, "chr")) == 0) {
  #   stop("It doesn't look deeptools matrix. Where are the genomic regions?")
  # }
  if (!startsWith(colnames(deeptools_matrix)[1], "@")) {
    stop("It doesn't look deeptools matrix.
      The first column name must start with @")
  }
  
  order = str_split(colnames(dt_mat)[1], "sample_labels")[[1]][2]
  print(glue('SAMPLE ORDER IN THE MATRIX: {order}'))
  
  # col 206: Mock center
  center_bins_mock = tibble("centerbin" = deeptools_matrix$V206,
                            type = "mock",
                            region = dt_mat$V4)
  # col 606: PDS center
  center_bins_pds = tibble("centerbin" = deeptools_matrix$V606,
                           type = "PDS",
                           region = dt_mat$V4)
  # col 1006: PhenDC3 center
  center_bins_phendc = tibble("centerbin" = deeptools_matrix$V1006,
                              type = "PhenDC3",
                              region = dt_mat$V4)
  center_bins = rbind(center_bins_mock, center_bins_pds, center_bins_phendc)
  
  if (plot) {
    p <- ggviolin(
      center_bins,
      x = "type",
      y = "centerbin",
      fill = "type",
      palette = mypal,
      add = "median_iqr"
    ) +
      labs(
        title = "",
        x = "",
        y = "centerbin",
        fill = " "
      ) +
      theme(
        text = element_text(size = 16),
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20, color = "black")
      ) + coord_cartesian(ylim = c(0, 10))
    annotate_figure(p, fig.lab = "", fig.lab.size = 6)
    print(p)
    
    ggsave(
      plot = p,
      filename = "../results/cutntag/centerbins-Sarma_CTCF_ChIP-seq-proms.pdf",
      device = "pdf",
      width = 4,
      height = 4
    )
    
  }
  
  return(center_bins)
}

# run
dt_mat = fread("../data/CTCF_all_with_promoters.matrix.gz")
extract_centerbin(dt_mat)

# save
extract_centerbin(dt_mat) %>% write_tsv("../results/cutntag/centerbins-Sarma_CTCF_ChIP-seq-proms.tsv")

df = extract_centerbin(dt_mat)
sc1 = df %>% dplyr::filter(type %in% c("mock", "PDS")) %>% 
  pivot_wider(names_from = type, values_from = centerbin)
sc1 = ggplot(sc1, aes(x = mock, y = PDS)) +
  geom_point(size = 2, color = "lightblue") +
  labs(
    title = "centerbin value - mock vs. PDS",
    x = "mock",
    y = "PDS",
  ) +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) 

sc2 = df %>% dplyr::filter(type %in% c("mock", "PhenDC3")) %>% 
  pivot_wider(names_from = type, values_from = centerbin)
sc2 = ggplot(sc2, aes(x = mock, y = PhenDC3)) +
  geom_point(size = 2, color = "lightblue") +
  labs(
    title = "centerbin value - mock vs. PhenDC3",
    x = "mock",
    y = "PhenDC3",
  ) +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_classic() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) 

scs = ggarrange(plotlist = list(sc1, sc2))
scs

ggsave(
  plot = scs,
  filename = "../results/cutntag/centerbins-Sarma_CTCF_ChIP-seq-proms_scatters.pdf",
  device = "pdf",
  width = 8,
  height = 4
)

# aggregate signals around the center bin
# helper function
around_centerbin = function(deeptools_matrix, centerbin, bin) {
  min_bin = centerbin - bin
  plus_bin = centerbin + bin
  mat = dt_mat[, min_bin:plus_bin]
  aggr = apply(mat, 1, max, na.rm=TRUE)
  return(aggr)
}

around_centerbin_30 = tibble(
  "mock" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 206,
    bin = 30
  ),
  "PDS" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 606,
    bin = 30
  ),
  "PhenDC3" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 1006,
    bin = 30
  )
)
around_centerbin_30 = pivot_longer(around_centerbin_30,
                                   "mock":"PhenDC3",
                                   names_to = "type",
                                   values_to = "aggr_around_centerbin")
around_centerbin_30 = around_centerbin_30 %>% mutate(bin_size = "30")

around_centerbin_40 = tibble(
  "mock" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 206,
    bin = 40
  ),
  "PDS" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 606,
    bin = 40
  ),
  "PhenDC3" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 1006,
    bin = 40
  )
)
around_centerbin_40 = pivot_longer(around_centerbin_40,
                                   "mock":"PhenDC3",
                                   names_to = "type",
                                   values_to = "aggr_around_centerbin")
around_centerbin_40 = around_centerbin_40 %>% mutate(bin_size = "40")

around_centerbin_60 = tibble(
  "mock" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 206,
    bin = 60
  ),
  "PDS" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 606,
    bin = 60
  ),
  "PhenDC3" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 1006,
    bin = 60
  )
)
around_centerbin_60 = pivot_longer(around_centerbin_60,
                                   "mock":"PhenDC3",
                                   names_to = "type",
                                   values_to = "aggr_around_centerbin")
around_centerbin_60 = around_centerbin_60 %>% mutate(bin_size = "60")

around_centerbin_80 = tibble(
  "mock" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 206,
    bin = 80
  ),
  "PDS" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 606,
    bin = 80
  ),
  "PhenDC3" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 1006,
    bin = 80
  )
)
around_centerbin_80 = pivot_longer(around_centerbin_80,
                                   "mock":"PhenDC3",
                                   names_to = "type",
                                   values_to = "aggr_around_centerbin")
around_centerbin_80 = around_centerbin_80 %>% mutate(bin_size = "80")


around_centerbin_100 = tibble(
  "mock" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 206,
    bin = 100
  ),
  "PDS" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 606,
    bin = 100
  ),
  "PhenDC3" = around_centerbin(
    deeptools_matrix = dt_mat,
    centerbin = 1006,
    bin = 100
  )
)
around_centerbin_100 = pivot_longer(around_centerbin_100,
                                    "mock":"PhenDC3",
                                    names_to = "type",
                                    values_to = "aggr_around_centerbin")
around_centerbin_100 = around_centerbin_100 %>% mutate(bin_size = "100")

around_centerbin_df = rbind(
  around_centerbin_30,
  around_centerbin_40,
  around_centerbin_60,
  around_centerbin_80,
  around_centerbin_100
)

p <- ggviolin(
  around_centerbin_df,
  x = "bin_size",
  y = "aggr_around_centerbin",
  fill = "type",
  palette = mypal,
  add = "median_iqr"
) +
  ylim(0, 10) +
  labs(
    title = "",
    x = "bin size around centerbin",
    y = "max aggregated signal",
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
#annotate_figure(p, fig.lab = "", fig.lab.size = 6)
print(p)

ggsave(
  plot = p,
  filename = "../results/cutntag/centerbins_around-Sarma_CTCF_ChIP-seq-proms.pdf",
  device = "pdf",
  width = 4,
  height = 4
)
