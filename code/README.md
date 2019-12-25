## Supporting code for Uncultivated Gut Virus manuscript

1. [Viral detection pipeline](viral_detection_pipeline/README.md): Identify viral sequences >=1Kb using the pipeline described in the manuscript

2. [Estimate genome completeness](genome_quality/README.md): Quantify genome completeness and apply genome standards

3. [Cluster genomes based on ANI](ani_cluster/README.md): Average nucleotide identity (ANI) code and centroid based clustering. Used to identify species-level clusters

4. [Cluster genomes based on AAI](ani_cluster/README.md). Average amino acid identity (AAI) code and MCL based clustering. Used to identify genus-level and family-level clusters

5. [Create marker-gene phylogenetic trees](marker_gene_tree/README.md). Identify prevalent single-copy genes in viral clade. Build phylogenetic tree based on concatenated alignments. Used to create phylogenies for family-level clusters

6. [Create SNP phylogenetic trees](snp_tree/README.md). Identify SNPs in core-genome regions based on whole-genome alignments. Build phylogenetic tree based on SNPs. Used to create strain-level phylogenies for species-level clusters


