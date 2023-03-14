#!/bin/bash
datafiles=$(ls Data/*.tsv | sed -e 's/Data\///' -e 's/\.tsv//')
for df in ${datafiles}; do 
    echo Calculating ssgsea pEMT for ${df}
    gseapy ssgsea -d Data/${df}.tsv -g ./Signatures/pEMT_gene_signature.gmt -o ./Output/${df}_pEMT --no-plot
done