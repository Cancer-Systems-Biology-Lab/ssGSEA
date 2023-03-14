#!/bin/bash
datafiles=$(ls Data/*.tsv | sed -e 's/Data\///' -e 's/\.tsv//')
for df in ${datafiles}; do 
    echo Calculating ssgsea EPI-MES for ${df}
    gseapy ssgsea -d Data/${df}.tsv -g ./Signatures/EM_gene_signature_cellLine_KS.gmt -o ./Output/${df}_EPIMES --no-plot
done