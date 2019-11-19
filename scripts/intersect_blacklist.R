require(tidyverse)
require(here)
require(glue)

# Which genes are within the blacklist?
bl_fname <- "hg38-blacklist.v2.bed"
bl_out <- here("nogit", bl_fname)
g_fname <- "umccr_cancer_genes.hg38.transcript.bed"
g_out <- here("nogit", g_fname)

readr::read_tsv(paste0("https://github.com/Boyle-Lab/Blacklist/raw/v2.0/lists/", bl_fname, ".gz"),
                col_names = c("chrom", "start", "end", "description"), col_types = "ciic") %>%
  dplyr::mutate(chrom = factor(chrom, levels = paste0("chr", c(1:22, "X", "Y")))) %>%
  dplyr::arrange(chrom, start) %>%
  dplyr::select(chrom, start, end) %>%
  readr::write_tsv(bl_out, col_names = FALSE)

readr::read_tsv(paste0("https://github.com/vladsaveliev/NGS_Utils/raw/433fe56586a2230a1059ad43e6e6888ae37ca5d8/",
                       "ngs_utils/reference_data/key_genes/", g_fname),
                col_names = c("chrom", "start", "end", "symbol"), col_types = "ciic") %>%
  dplyr::mutate(chrom = factor(chrom, levels = paste0("chr", c(1:22, "X", "Y")))) %>%
  dplyr::arrange(chrom, start) %>%
  readr::write_tsv(g_out, col_names = FALSE)

bedtools <- "/Users/pdiakumis/my_apps/miniconda/envs/woof/bin/bedtools"
query <- glue::glue("{bedtools} intersect -a {bl_out} -b {g_out} -wb")
system(query, intern = TRUE) %>%
  tibble::tibble(all_cols = .) %>%
  tidyr::separate(col = .data$all_cols,
                  into = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "symbol"),
                  sep = "\t", convert = TRUE)
# # A tibble: 17 x 7
# chrom1    start1      end1 chrom2    start2      end2 symbol
# <chr>      <int>     <int> <chr>      <int>     <int> <chr>
# 1 chr3    36993331  36995000 chr3    36993331  37050918 MLH1
# 2 chr3    49683946  49689053 chr3    49683946  49689053 MST1
# 3 chr3    78948800  78950500 chr3    78597239  79767815 ROBO1
# 4 chr3    89345600  89370500 chr3    89107523  89482134 EPHA3
# 5 chr3   195746764 195750100 chr3   195746764 195812277 MUC4
# 6 chr3   195775500 195791400 chr3   195746764 195812277 MUC4
# 7 chr6      292096    351355 chr6      292096    351355 DUSP22
# 8 chr7   152375800 152435100 chr7   152134928 152436005 KMT2C
# 9 chr9      319900    322400 chr9      214864    465259 DOCK8
# 10 chr12     371800    389454 chr12     280128    389454 KDM5A
# 11 chr12  113079600 113081500 chr12  113057689 113098028 DTX1
# 12 chr12  124430500 124440300 chr12  124324414 124535603 NCOR2
# 13 chr18   47852300  47854300 chr18   47831550  47931146 SMAD2
# 14 chr18   52791800  52793800 chr18   52340171  53535899 DCC
# 15 chrX     1200600   1209400 chrX     1190489   1212750 CRLF2
# 16 chrX     1480900   1492800 chrX     1462571   1537107 P2RY8
# 17 chrX    67626800  67632300 chrX    67544035  67730619 AR
