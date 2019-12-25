# SNP-based phylogenetic tree

<b> Requirements </b>  

* MUMmer v4.0.0beta2
* FastTreeMP v2.1.10

You should be able to run the following from your command line:
`nucmer`  
`FastTreeMP`  

<b> Download viral genomes from species-cluster 51243 </b>  
This is the Faecalibacterium prauznitzi phage highlighted in Fig 3. Other input data can be used  
`wget https://www.dropbox.com/s/lblvqunl9lmsbpk/OTU-51243.tar.gz && tar -zxvf OTU-51243.tar.gz`  

<b> Align all genomes to a selected representative </b>  
`python align_genomes.py --genomes OTU-51243 --ref OTU-51243/UGV-GENOME-0380053.fna --out alignments`  

<b> Use SNPs to build multiple sequence alignment </b>  
`python build_msa.py --in alignments --out snps.fna --max_gaps_col 50 --max_gaps_seq 50`  

<b> Use FastTree to contruct phylogeny </b>  
`FastTreeMP snps.fna > snps.tree`  

<b> Pipeline overview </b>  

1. Use `nucmer` utility in MUMmer4 to align all genomes to a reference
2. Call SNPs in 1:1 alignment blocks using `show-snps, show-coords, and show-diff` utlilties
3. Use alignments to create multiple sequence alignments against reference
4. Remove genomic sites covered in <50% of genomes and non-polymorphic sites (i.e. non-SNPs)
5. Remove genomes containing >50% gaps
6. Use `FastTree` to create phylogeny from trimmed multiple sequence alignment


