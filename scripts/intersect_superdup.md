---
title: "Overlapping Genes with SuperDups"
author: "Peter Diakumis"
date: "Tue 2019-Nov-26"
output:
  html_document:
    keep_md: yes
---




```r
require(tidyverse)
require(here)
require(glue)
```

## bedtools intersect


```r

tx_hg19 <- here("nogit/transcripts/ngs_utils/umccr_cancer_genes.GRCh37.transcript.bed")
tx_hg38 <- here("nogit/transcripts/ngs_utils/umccr_cancer_genes.hg38.transcript.bed")
cod_hg19 <- here("nogit/transcripts/ngs_utils/umccr_cancer_genes.GRCh37.coding.bed")
cod_hg38 <- here("nogit/transcripts/ngs_utils/umccr_cancer_genes.hg38.coding.bed")

# cut -f1-3 GRCh37GenomicSuperDup.tab | grep -v "gl" | grep -v "chrom" | sort -k1,1V -k2,2n -k3,3n | sed 's/chr//' | uniq > GRCh37GenomicSuperDup_clean.bed
# cut -f1-3 GRCh38GenomicSuperDup.tab | sort -k1,1V -k2,2n -k3,3n | grep -v "_" | grep -v "EBV" | uniq > GRCh38GenomicSuperDup_clean.bed
dup_hg19 <- here("nogit/superdups/GRCh37GenomicSuperDup_clean.bed")
dup_hg38 <- here("nogit/superdups/GRCh38GenomicSuperDup_clean.bed")

bedtools_intersect <- function(a, b) {
  bedtools <- "/Users/pdiakumis/my_apps/miniconda/envs/woof/bin/bedtools"
  query <- glue::glue("{bedtools} intersect -a {a} -b {b} -wo")
  system(query, intern = TRUE) %>%
    tibble::tibble(all_cols = .) %>%
    tidyr::separate(col = .data$all_cols,
                    into = c("chrom1", "start1", "end1", "symbol", "chrom2", "start2", "end2", "overlap"),
                    sep = "\t", convert = TRUE)
}

(hg19_tx_dup <- bedtools_intersect(tx_hg19, dup_hg19))
## # A tibble: 690 x 8
##    chrom1   start1     end1 symbol chrom2   start2     end2 overlap
##    <chr>     <int>    <int> <chr>  <chr>     <int>    <int>   <int>
##  1 1       6245079  6259672 RPL22  1       6245088  6246878    1790
##  2 1       6245079  6259672 RPL22  1       6245089  6246878    1789
##  3 1       9711802  9788977 PIK3CD 1       9787004  9788973    1969
##  4 1      15817326 15850940 CASP9  1      15799949 15818262     936
##  5 1      17733255 17766220 RCC2   1      17733257 17735335    2078
##  6 1      17733255 17766220 RCC2   1      17733257 17735691    2434
##  7 1      17733255 17766220 RCC2   1      17733258 17735692    2434
##  8 1      17733255 17766220 RCC2   1      17733275 17735692    2417
##  9 1      17733255 17766220 RCC2   1      17733276 17735693    2417
## 10 1      22379119 22419437 CDC42  1      22417917 22419436    1519
## # … with 680 more rows
(hg38_tx_dup <- bedtools_intersect(tx_hg38, dup_hg38))
## # A tibble: 722 x 8
##    chrom1   start1     end1 symbol chrom2   start2     end2 overlap
##    <chr>     <int>    <int> <chr>  <chr>     <int>    <int>   <int>
##  1 chr1    6185019  6199612 RPL22  chr1    6185028  6186818    1790
##  2 chr1    6185019  6199612 RPL22  chr1    6185029  6186818    1789
##  3 chr1    9651744  9728919 PIK3CD chr1    9726946  9728915    1969
##  4 chr1   15490831 15524445 CASP9  chr1   15473454 15491767     936
##  5 chr1   17408675 17438561 RCC2   chr1   17406761 17408839     164
##  6 chr1   17408675 17438561 RCC2   chr1   17406761 17409195     520
##  7 chr1   17408675 17438561 RCC2   chr1   17406762 17409196     521
##  8 chr1   17408675 17438561 RCC2   chr1   17406779 17409196     521
##  9 chr1   17408675 17438561 RCC2   chr1   17406780 17409197     522
## 10 chr1   22052626 22092946 CDC42  chr1   22086436 22092405    5969
## # … with 712 more rows
(hg19_cod_dup <- bedtools_intersect(cod_hg19, dup_hg19))
## # A tibble: 1,683 x 8
##    chrom1   start1     end1 symbol chrom2   start2     end2 overlap
##    <chr>     <int>    <int> <chr>  <chr>     <int>    <int>   <int>
##  1 1       6246731  6246876 RPL22  1       6245088  6246878     145
##  2 1       6246731  6246876 RPL22  1       6245089  6246878     145
##  3 1       9786966  9787104 PIK3CD 1       9787004  9788973     100
##  4 1      17735585 17735690 RCC2   1      17733257 17735691     105
##  5 1      17735585 17735690 RCC2   1      17733258 17735692     105
##  6 1      17735585 17735690 RCC2   1      17733275 17735692     105
##  7 1      17735585 17735690 RCC2   1      17733276 17735693     105
##  8 1      22412931 22413041 CDC42  1      22412929 22418898     110
##  9 1      22412931 22413041 CDC42  1      22412929 22418907     110
## 10 1      22412931 22413041 CDC42  1      22412929 22419443     110
## # … with 1,673 more rows
(hg38_cod_dup <- bedtools_intersect(cod_hg38, dup_hg38))
## # A tibble: 1,794 x 8
##    chrom1   start1     end1 symbol chrom2   start2     end2 overlap
##    <chr>     <int>    <int> <chr>  <chr>     <int>    <int>   <int>
##  1 chr1    6186671  6186816 RPL22  chr1    6185028  6186818     145
##  2 chr1    6186671  6186816 RPL22  chr1    6185029  6186818     145
##  3 chr1    9726908  9727046 PIK3CD chr1    9726946  9728915     100
##  4 chr1   17409089 17409194 RCC2   chr1   17406761 17409195     105
##  5 chr1   17409089 17409194 RCC2   chr1   17406762 17409196     105
##  6 chr1   17409089 17409194 RCC2   chr1   17406779 17409196     105
##  7 chr1   17409089 17409194 RCC2   chr1   17406780 17409197     105
##  8 chr1   22086438 22086548 CDC42  chr1   22086436 22092405     110
##  9 chr1   22086438 22086548 CDC42  chr1   22086436 22092414     110
## 10 chr1   22086438 22086548 CDC42  chr1   22086436 22092950     110
## # … with 1,784 more rows
```

### Overlap summary


```r
purrr::map_dfr(list(hg19_tx_dup = hg19_tx_dup, hg38_tx_dup = hg38_tx_dup,
                    hg19_cod_dup = hg19_cod_dup, hg38_cod_dup = hg38_cod_dup),
           function(d) {
             pull(d, overlap) %>% summary() %>% as.matrix() %>% t() %>% as_tibble()
       }, .id = "pair")
## # A tibble: 4 x 7
##   pair          Min. `1st Qu.` Median  Mean `3rd Qu.`   Max.
##   <chr>        <dbl>     <dbl>  <dbl> <dbl>     <dbl>  <dbl>
## 1 hg19_tx_dup    196     1398.  2476. 7193.     7134.  88139
## 2 hg38_tx_dup    158     1433   2508  7994.     7390  169199
## 3 hg19_cod_dup     3       95    129   202.      174    5793
## 4 hg38_cod_dup     1       95    131   200.      176    5793
```

### Number of genes


```r
purrr::map_int(list(hg19_tx_dup = hg19_tx_dup, hg38_tx_dup = hg38_tx_dup,
                    hg19_cod_dup = hg19_cod_dup, hg38_cod_dup = hg38_cod_dup),
               function(d) {
                 pull(d, symbol) %>% table() %>% length()
       })
##  hg19_tx_dup  hg38_tx_dup hg19_cod_dup hg38_cod_dup 
##          221          219          139          139
```

