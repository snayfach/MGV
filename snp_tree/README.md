# SNP-based phylogenetic tree

<b> Requirements </b>  

* MUMmer v4.0.0beta2
* FastTreeMP v2.1.10

You should be able to run the following from your command line:  
`nucmer`  
`FastTreeMP`  

<b> Locate viral genomes for an individual viral OTU and place these into a directory. Choose one of these genomes to be the representative </b>     
GENOMES=`my_genomes`  
REP=`my_genomes/representative.fna`  

<b> Align all genomes to a selected representative </b>  
`python align_genomes.py --genomes $GENOMES --ref $REP --out alignments`  

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


