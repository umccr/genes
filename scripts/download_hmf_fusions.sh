#!/usr/bin/env bash

# links work as of 2019-Oct-02

base_url="https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMFTools-Resources%2FLINX&files="
wget ${base_url}knownPromiscuousFive.csv -O hmf_knownPromiscuousFive.csv
wget ${base_url}knownPromiscuousThree.csv -O hmf_knownPromiscuousThree.csv
wget ${base_url}knownFusionPairs.csv -O hmf_knownFusionPairs.csv
