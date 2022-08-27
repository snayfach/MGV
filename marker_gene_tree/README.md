# Marker gene phylogenetic tree

Note: this code was not used to create the Caudovirales phylogeny shown in figure 4. For that, please follow the description of methods from https://doi.org/10.1038/s41564-019-0448-z. This code is useful to create a phylogeny for an order or family level viral clade. But it is not sensitive enough to find sufficient prevalent genes across an entire viral class like the Caudoviricites

<b> Requirements </b>  

* HMMER v3.1b2
* FAMSA v1.2.5
* trimAL v1.4.rev22
* FastTreeMP v2.1.10
* MCL v14-137

You should be able to run the following from your command line:  
`hmmsearch`  
`famsa`  
`trimal`  
`FastTreeMP`  
`mcl`  
`prodigal`  

<b> Locate viral proteins of genomes from an individual clade </b>  

* proteins should have headers with the format: `[GENOME_ID]_[PROTEIN_NUM]`
* this is the format that Prodigal will produce  

<b> Run pipeline </b>  
`python marker_tree.py --in_faa proteins.faa --out_dir out`

<b> Pipeline overview </b>  

1. Make diamond db of proteins
2. Perform all-vs-all alignment of proteins using diamond (E-value<1e-5 )
3. Run MCL to cluster proteins (1.4 inflation factor)
4. Select protein clusters (copy number <1.1, prevalence >10%)
5. Write seqs for each PC
6. Build multiple seq alignment (MSA) for each PC using FAMSA
7. Build hmm from each MSA using HMMER
8. Search hmms versus original proteins (E-value<1e-5 )
9. Extract proteins for top hmm hits
10. Build multiple seq alignment (MSA) for each PC using expanded set of homologs using FAMSA
11. Trim multiple seq alignments with trimAl (discard columns with >50% gaps)
12. Concatenate alignments and fill missing genes with gaps
13. Count gaps per genome
14. Write concatenated alignment for genomes with <90% gaps
15. Built phylogeny with FastTreeMP using default options


