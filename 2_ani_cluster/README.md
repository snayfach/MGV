# Cluster viral genomes based on ANI

<b> Locate viral genomes to cluster </b>
FNA=../0_viral_detect/output/SRS1735492.fna

<b> Perform all-vs-all BLAST </b>
`makeblastdb -in $FNA -out blastdb -dbtype nucl`
`blastn -query $FNA -db blastdb -out blast.tsv -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90`

<b> 2. Compute ANI from BLAST results </b>
`python blastani.py -i blast.tsv -o ani.tsv`

<b> Perform centroid-based clustering </b>
`python cluster.py --fna $FNA --ani ani.tsv --out clusters.tsv --min_ani 95 --min_qcov 0 --min_tcov 50`

