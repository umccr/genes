#!/usr/bin/env bash

# link works as of 2019-Oct-02
# seds by Vlad Saveliev

wget https://raw.githubusercontent.com/ndaniel/fusioncatcher/master/bin/generate_known.py -O - | \
  grep "        \['" | \
  sed "s#        \['##" | \
  sed "s#','#,#" | \
  sed "s#'\],##" | \
  sed "s#'\]##" | \
  sort -u > fusioncatcher_pairs.txt