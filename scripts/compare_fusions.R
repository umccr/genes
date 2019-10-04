# Author: Vlad Saveliev
# Maintainer: Peter Diakumis
################ Fusions #################
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(purrr)


# simple_sv_annotation list is coming from FusionCatcher,
# so we get it from there to make sure it's most recent
# (e.g. 11-Feb-2019: 7866, 02-Oct-2019: 8390 pairs versus 6527 in simple_sv_annotation):
# H_gene = Head = 5'
# T_gene = Tail = 3'

# 8,391 pairs
fus_catcher <- readr::read_tsv("../fusions/fusioncatcher_pairs.txt",
                        col_names = "pair", col_types = "c") %>%
  tidyr::separate(pair, c("H_gene", "T_gene"), sep = ",") %>%
  dplyr::distinct()


# 393 pairs
hmf_pairs <-
  readr::read_csv("../fusions/hmf_knownFusionPairs.csv", col_types = "cc", col_names = TRUE) %>%
  dplyr::rename(H_gene = fiveGene, T_gene = threeGene)

hmf_prom_head <- read_csv("../fusions/hmf_knownPromiscuousFive.csv", col_types = "c", col_names = TRUE) # 28 genes
hmf_prom_tail <- read_csv("../fusions/hmf_knownPromiscuousThree.csv", col_types = "c", col_names = TRUE) # 37 genes

# 398 fusion genes
cancer_genes <- read_tsv("../key_genes/umccr_cancer_genes.2019-03-20.tsv") %>%
  filter(fusion == T) %>%
  select(gene = symbol)

pairs <- hmf_pairs %>%
  full_join(fus_catcher %>% mutate(fus_catch = TRUE), by = c("H_gene", "T_gene"))

# Are there genes that are both head and tail?
pairs %>% filter(H_gene %in% pairs$T_gene) %>% distinct(H_gene)
# 1318 such genes, e.g.:
pairs %>% filter(H_gene == 'ACPP' | T_gene == 'ACPP')
# and tend to be in the fusioncatcher list, which is gigantic compared to HMF. Try to remove from it promiscuous fusions?


# Do HMF pairs cover fusioncatcher?
fus_catcher %>%
  unite(fus, H_gene, T_gene, sep = '&') %>%
  filter(fus %in% (unite(hmf_pairs, fus, H_gene, T_gene, sep = '&')$fus))
# 255 / 7632 (down to 227 / 8391 now)
fus_catcher %>%
  unite(fus, H_gene, T_gene, sep = '&') %>%
  filter(fus %in% (unite(hmf_pairs, fus, T_gene, H_gene, sep = '&')$fus))
# - plus 174 (200 now) if we swap T and H
# Do fusioncatcher cover HMF pairs?
hmf_pairs %>%
  unite(fus, H_gene, T_gene, sep = '&') %>%
  filter(fus %in% (unite(fus_catcher, fus, H_gene, T_gene, sep = '&')$fus))
# 255 / 401 (down to 227 / 398 now)

# Do HMF promiscous cover fusioncatcher?
fus_catcher %>%
  filter(H_gene %in% hmf_prom_head$gene | T_gene %in% hmf_prom_tail$gene | T_gene %in% hmf_prom_head$gene | H_gene %in% hmf_prom_tail$gene)
# 1659 / 7632 (1953 / 8391)
# So promiscuous cover only 14% of fusions, so we better stick to HMF fusions only.
# Also, do fusioncatcher cover HMF promiscous?
hmf_prom_head %>% filter(gene %in% fus_catcher$H_gene | gene %in% fus_catcher$T_gene)
# 29/30 (28 / 28 now)
hmf_prom_tail %>% filter(gene %in% fus_catcher$H_gene | gene %in% fus_catcher$T_gene)
# 36/36 (37 / 37 now)

# https://github.com/pmelsted/pizzly/issues/19
fus_catcher %>% filter(str_detect(H_gene, "IGH\\.*"))
fus_catcher %>% filter(str_detect(H_gene, "DUX4"))

################
# How about cancer genes?
cancer_genes
# 352 (398 now)
hmf_pairs %>%
  count(H_gene %in% cancer_genes$gene, T_gene %in% cancer_genes$gene)
#`H_gene %in% cancer_genes$gene` `T_gene %in% cancer_genes$gene`     n
# FALSE                           FALSE                              49
# FALSE                           TRUE                               72
# TRUE                            FALSE                              62
# TRUE                            TRUE                              218
# 218 full pairs, 62 heads only, 72 tails only match, and only 49 pairs do not match completely.
# Now:
#`H_gene %in% cancer_genes$gene` `T_gene %in% cancer_genes$gene`     n
#1 FALSE                           FALSE                            29
#2 FALSE                           TRUE                             66
#3 TRUE                            FALSE                            51
#4 TRUE                            TRUE                            247

# Also:
hmf_prom_head %>% count(gene %in% cancer_genes$gene)  # 25/36 (28/28)
hmf_prom_tail %>% count(gene %in% cancer_genes$gene)  # 27/30 (35/37)
# Mostly matching, so we'll just add remaining fusions in the cancer gene list, and use HMF list of fusions in simple_sv_annotation.
