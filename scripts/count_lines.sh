#!/usr/bin/env bash

printf "#### Fusion counts\n\n"
printf "\`\`\`\n"

for f in $(find ../fusions -type f)
do
    l=$(wc -l $f | awk '{print $1'})
    printf "$f: $l\n"
done

printf "\`\`\`\n"