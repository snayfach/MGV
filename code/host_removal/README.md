# Remove host regions

<b> Script to remove host regions </b>  
`remove_host.py --in_base <BASENAME> --out_base <BASENAME>`

Expected input files (output files from viral detection pipeline):  

* BASENAME.tsv # viral features table
* BASENAME.fna
* BASENAME.ffn
* BASENAME.faa
* BASENAME.gff

The same set of files will be written to --out_base


<b> Pipeline overview </b>  

1. Scan each viral contig, starting from the 5' and 3' ends
2. Identify the largest block of genes, containing only Pfam-annotated or unannotated genes, but no viral-annotated genes, and ending with a Pfam-annotated gene
3. If the gene block contains at least 2 Pfam-annotated genes and at least 30% of genes have a Pfam-annotation, then the gene block is determined to be a host region
4. The genomic coordinates of the block are determined, and removed from input sequences


