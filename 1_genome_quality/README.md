# Estimate genome quality and completeness

<b>  Necessary dependencies </b>  
DIAMOND   

You should be able to run the following from the command line:  
`diamond`  

<b> Move to the correct subdirectory in repository</b>  
`cd /path/to/UGV/1_viral_qc`  

<b> Locate input viral sequences </b>
FNA=../0_viral_detect/output/SRS1735492.fna
FAA=../0_viral_detect/output/SRS1735492.faa
FFN=../0_viral_detect/output/SRS1735492.ffn

<b> Download database of complete or circular reference genomes </b>
* Data sources: IMG/VR (16556 circular), UGV (30591 circular), NCBI GenBank (30329 complete)
`wget https://www.dropbox.com/s/8af4g2k4ec011a2/complete_refs_NCBI_IMG_UGV.tar.gz && tar -zxvf complete_refs_NCBI_IMG_UGV.tar.gz`

<b> Estimate completeness of viral contigs </b>
* Find nearest reference with >80% AAI over >50% genes to each contig
* Completeness based on ratio of sequence lengths
`python completeness.py --in_fna $FNA --in_faa $FAA --out_m8 blastp.tsv --out_tsv completeness.tsv --threads 10 --min_gene_sharing 50 --min_aai 80 --db_dir complete_refs_NCBI_IMG_UGV`

<b> Identify circular contigs </b>
* based on a >=20 bp non-inverted terminal repeat
`python circularity.py --in_fna $FNA --out_tsv circularity.tsv --repeat_length 20`
