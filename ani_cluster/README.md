# Cluster viral genomes based on ANI

<b> Locate viral genomes to cluster </b>  
In the README, we've specified the nucleotide sequences from viruses identified from SRS1735492, but any sequences can be used

<b> Perform all-vs-all BLAST using megablast utility </b>  
`makeblastdb -in SRS1735492.fna -out blastdb -dbtype nucl`
`blastn -query SRS1735492.fna -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90`

Note that the default utility in blastn is megablast. We indicate that we wish to report many alignments (up to 25000) to make sure that we estimate ANI between all genome pairs. We also indicate that we only want alignments above 90% identity since we are only looking for viruses that are within 95% identity. This reduces the size of the output quite a bit.

<b> 2. Compute ANI from BLAST results </b>  
`python blastani.py -i blast.tsv -o ani.tsv`

Average nucleotide identity is based on the length-weighted average DNA identity across all local alignments between each pair of genomes. This ANI estimate is comparable to MUMmer and MiSi, but much faster. The alignment fraction (AF) is computed based on the length of merged alignment coordinates relative to each genome.

<b> Perform centroid-based clustering </b>  
`python cluster.py --fna SRS1735492.fna --ani ani.tsv --out clusters.tsv --min_ani 95 --min_qcov 0 --min_tcov 85`

Clustering is performed using a greedy, centroid-based algorithm in which:

1. sequences are sorted by length  
2. the longest contig is designated as the centroid of a new cluster
3. all contigs within 95% ANI and 85% AF are assigned to that cluster
4. steps (2-3) are repeated until all sequences had been assigned to a cluster
