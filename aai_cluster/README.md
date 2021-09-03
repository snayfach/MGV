# Cluster genomes based on AAI

<b> Locate viral proteins from genomes</b>  
In the README, we've specified the path to proteins from viruses identified from SRS1735492, but any sequences can be used. The code assumes each sequence header has the format '>CONTIG-ID_GENE-NUM'   

<b> Make DIAMOND database </b>  
`diamond makedb --in SRS1735492.faa --db viral_proteins --threads 10`

<b> Perform all-vs-all BLASTP </b>  
`diamond blastp --query SRS1735492.faa --db viral_proteins --out blastp.tsv --outfmt 6 --evalue 1e-5 --max-target-seqs 10000 --query-cover 50 --subject-cover 50`

<b> Compute AAI from BLAST results </b>   
`python amino_acid_identity.py --in_faa SRS1735492.faa --in_blast blastp.tsv --out_tsv aai.tsv`  

Amino acid identity is computed based on the average BLAST percent identity between all genes shared between each pair of genomes (E-value <1e-5)  

<b> Filter edges and prepare MCL input </b>    
`python filter_aai.py --in_aai aai.tsv --min_percent_shared 20 --min_num_shared 16 --min_aai 40 --out_tsv genus_edges.tsv`  
`python filter_aai.py --in_aai aai.tsv --min_percent_shared 10 --min_num_shared 8 --min_aai 20 --out_tsv family_edges.tsv`  

Here we're keeping edges between genomes with >=20% AAI and genomes with either 8 shared genes or at least 20% of shared genes (relative to both genomes)

<b> Perform MCL-based clustering </b>  
`mcl genus_edges.tsv -te 8 -I 2.0 --abc -o genus_clusters.txt`
`mcl family_edges.tsv -te 8 -I 1.2 --abc -o family_clusters.txt`

In the output each row indictes the members belonging to each cluster (including singletons)

