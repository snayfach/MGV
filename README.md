## Supporting materials for manuscript: "Metagenomic compendium of 189,680 DNA viruses from the human gut microbiome"

### Data availability

See: https://portal.nersc.gov/MGV/

### Code availability

1. [Viral detection pipeline](viral_detection_pipeline/README.md): Identify viral contigs >=1Kb using the pipeline described in the manuscript

2. [Quality control](https://bitbucket.org/berkeleylab/checkv): Identify and remove putative host regions flanking viral contigs. Quantify genome completeness and apply genome quality standards.

3. [Cluster genomes based on ANI](ani_cluster/README.md): Average nucleotide identity (ANI) code and centroid based clustering. Used to identify species-level viral clusters

4. [Cluster genomes based on AAI](aai_cluster/README.md). Average amino acid identity (AAI) code and MCL based clustering. Used to identify genus-level and family-level viral clusters

5. [Create SNP phylogenetic trees](snp_tree/README.md). Identify SNPs in core-genome regions based on whole-genome alignments. Build phylogenetic tree based on SNPs. Used in manuscript to create strain-level phylogenies for species-level viral clusters

6. [Create marker-gene phylogenetic trees](marker_gene_tree/README.md). Identify prevalent single-copy genes in a viral clade. Use concatenated gene alignments to build phylogenetic tree. 

7. [Identify CRISPR spacers](crispr_spacers/README.md). Identify CRISPR spacers using CRT and PILERCR, merge redundant CRISPR arrays, and format output.

For any other code/analysis inquiries, please open a github issue. Note: most of these scripts were written for Python 2. If you get an error using Python 3, try re-running with Python 2.

If this code is useful, please cite:
Nayfach et al. Metagenomic compendium of 189,680 DNA viruses from the human gut microbiome. 2021. https://www.nature.com/articles/s41564-021-00928-6.  

### Updated version of the database
See: https://github.com/snayfach/UHGV



