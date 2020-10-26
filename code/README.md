## Supporting code for Uncultivated Gut Virus manuscript

1. [Viral detection pipeline](viral_detection_pipeline/README.md): Identify viral contigs >=1Kb using the pipeline described in the manuscript

2. [Quality control](https://bitbucket.org/berkeleylab/checkv): Identify and remove putative host regions flanking viral contigs. Quantify genome completeness and apply genome quality standards.

3. [Cluster genomes based on ANI](ani_cluster/README.md): Average nucleotide identity (ANI) code and centroid based clustering. Used to identify species-level viral clusters

4. [Cluster genomes based on AAI](ani_cluster/README.md). Average amino acid identity (AAI) code and MCL based clustering. Used to identify genus-level and family-level viral clusters

5. [Create marker-gene phylogenetic trees](marker_gene_tree/README.md). Identify prevalent single-copy genes in a viral clade. Use concatenated gene alignments to build phylogenetic tree. Used in manuscript to create phylogenies for family-level viral clusters

6. [Create SNP phylogenetic trees](snp_tree/README.md). Identify SNPs in core-genome regions based on whole-genome alignments. Build phylogenetic tree based on SNPs. Used in manuscript to create strain-level phylogenies for species-level viral clusters

7. [Identify CRISPR spacers](crispr_spacers/README.md). Identify CRISPR spacers using CRT and PILERCR, merge redundant CRISPR arrays, and format output.

