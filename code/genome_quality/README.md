# Estimate genome quality and completeness

<b>  Necessary dependencies </b>  
DIAMOND   

You should be able to run the following from the command line:  
`diamond`  

<b> Move to the correct subdirectory in repository</b>  
`cd /path/to/UGV/1_viral_qc`  

<b> Locate input viral sequences </b>  
In the README, we've specified the sequences from viruses identified from SRS1735492, but any sequences can be used

<b> Download database of complete or circular reference genomes </b>  
`wget https://www.dropbox.com/s/8af4g2k4ec011a2/complete_refs_NCBI_IMG_UGV.tar.gz && tar -zxvf complete_refs_NCBI_IMG_UGV.tar.gz`  

These sequences were derived from the following data sources: IMG/VR (16,556 circular), UGV (30,591 circular), NCBI GenBank (30,329 complete). Circular sequences are at least 2Kbp and contain a 20-bp terminal repeat.

<b> Estimate completeness of viral contigs </b>  
`python completeness.py --in_fna SRS1735492.fna --in_faa SRS1735492.faa --out_m8 blastp.tsv --out_tsv completeness.tsv --threads 10 --min_gene_sharing 50 --min_aai 80 --db_dir complete_refs_NCBI_IMG_UGV`

This code works by:  
1. Aligning proteins from query genomes versus references  
2. Computing the average amino acid identity to references (AAI)  
3. For each query, identify the nearest reference that has >80% AAI and shares >50% of genes with query  
4. Computing the ratio of genome lengths: 100.0*query/reference  
5. A value of 50% indicate the query is 1/2 the expected length  
6. Maximum reported values are 99%  

<b> Identify circular contigs </b>
`python circularity.py --in_fna SRS1735492.fna --out_tsv circularity.tsv --repeat_length 20`

Circularity is determined based on a >=20 bp non-inverted terminal repeat.

<b> Criterea for selecting draft-quality genomes </b>  
1. Circular and >=2Kbp, OR  
2. >=50% complete and >=2Kbp, OR  
3. >=10Kbp where non-circular and completeness is not known  