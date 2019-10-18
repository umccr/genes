Genes
=====
Collections of genes used at UMCCR.
- [Genes](#genes)
  - [Known fusions](#known-fusions)

## CPSR Cancer Predisposition Genes
Includes 209 genes from [CPSR](https://github.com/sigven/cpsr/blob/master/predisposition.md),
plus a further 8 genes based on advice from curators. Total: 217 genes (2019-Oct-10).











## Known fusions
Known fusions are curated by HMF. From their
[2018 preprint](https://www.biorxiv.org/content/biorxiv/early/2018/09/20/415133.full.pdf):

```
6. Identification of gene fusions

For each structural variant, every combination of annotated overlapping
transcripts from each breakend was tested to see if it could potentially
form an intronic inframe fusion. A list of 411 curated known fusion pairs
was sourced by taking the union of known fusions from the following
external databases:

- Cosmic curated fusions (v83)
- OncoKb (download = 01-mar-2018)
- CGI (update: 17-jan-2018)
- CIViC (download = 01-mar-2018)

We then also created a list of promiscuous fusion partners using the
following rules:

- 3'/5' promiscuous: Any gene which appears on the [3'/5'] side in more than 3
of the curated fusion pairs OR appears at least once on the [3'/5'] side and
is marked as promiscuous in either OncoKb, CGI or CIVIC

For each promiscuous partner we also curated a list of essential
domains that must be preserved to form a viable fusion partner.
Finally, we report an intronic inframe fusion if the following conditions
are met:

- matches an exact fusion from the curated list OR
- is intergenic and matches 5’ promiscuous OR
- matches 3’ promiscuous gene

```

To download these fusions (as of 2019-Oct-02):

```
base_url="https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMFTools-Resources%2FLINX&files="
wget ${base_url}knownPromiscuousFive.csv -O hmf_knownPromiscuousFive.csv
wget ${base_url}knownPromiscuousThree.csv -O hmf_knownPromiscuousThree.csv
wget ${base_url}knownFusionPairs.csv -O hmf_knownFusionPairs.csv
```

There is also a more broad list of fusions from [FusionCatcher](https://github.com/ndaniel/fusioncatcher).
This list is within a Python script. Yes.

To download these fusions (as of 2019-Oct-02):

```
wget https://raw.githubusercontent.com/ndaniel/fusioncatcher/master/bin/generate_known.py -O - | \
  grep "        \['" | \
  sed "s#        \['##" | \
  sed "s#','#,#" | \
  sed "s#'\],##" | \
  sed "s#'\]##" | \
  sort -u > fusioncatcher_pairs.txt
```

There are also known fusion genes in [NGC](http://ncg.kcl.ac.uk/), one of the sources for UMCCR cancer gene list.
We compare the lists with `compare.R`.
The FusionsCatcher list is too big so we assign a lower tier to a matching fusion when prioritizing.
