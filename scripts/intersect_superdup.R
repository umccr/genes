require(tidyverse)
require(here)
require(glue)

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
  query <- glue::glue("{bedtools} intersect -a {a} -b {b} -wb")
  system(query, intern = TRUE) %>%
    tibble::tibble(all_cols = .) %>%
    tidyr::separate(col = .data$all_cols,
                    into = c("chrom1", "start1", "end1", "symbol", "chrom2", "start2", "end2"),
                    sep = "\t", convert = TRUE)
}

(hg19_tx_dup <- bedtools_intersect(tx_hg19, dup_hg19))
(hg38_tx_dup <- bedtools_intersect(tx_hg38, dup_hg38))
(hg19_cod_dup <- bedtools_intersect(cod_hg19, dup_hg19))
(hg38_cod_dup <- bedtools_intersect(cod_hg38, dup_hg38))
