require(tidyverse)
require(rvest)
require(xml2)
require(here)

# Date: 2019-10-18, CPSR v0.5.1, commit hash shown in URL
# Total of 209 genes
cpsr_url <- "https://raw.githubusercontent.com/sigven/cpsr/96613b91b3fec194fa4b8d1249f8d631bab2e4e2/predisposition_genes_20181112.tsv"
cpsr_genes_2019_10_18 <- readr::read_tsv(cpsr_url, col_names = T) %>%
  dplyr::distinct() %>%
  dplyr::arrange(symbol)

# Total of 216 genes
predispose_genes <-
  readr::read_tsv(here("key_genes/sources/pcgr_predispose_genes.txt"), col_names = "symbol", col_types = "c") %>%
  dplyr::distinct()

table(cpsr_genes_2019_10_18$symbol %in% predispose_genes$symbol) # 1, 208
table(predispose_genes$symbol %in% cpsr_genes_2019_10_18$symbol) # 8, 208

# which genes are in predispose list, but not CPSR one?
predispose_genes[!predispose_genes$symbol %in% cpsr_genes_2019_10_18$symbol, ]
# c("TNFRSF6", "KLLN", "MAP3K6", "NEK1", "NTRK1", "RAD54L", "RHNO1", "RTEL1")

# which genes are in CPSR list, but not predispose one?
cpsr_genes_2019_10_18$symbol[!cpsr_genes_2019_10_18$symbol %in% predispose_genes$symbol]
# c("FAS")

# include all 217 (216 + 1) into one
genes_all <- unique(c(predispose_genes$symbol, cpsr_genes_2019_10_18$symbol)) %>% sort()

readr::write_lines(genes_all, here("sources/cpsr_cancer_predisposition.tsv"))
