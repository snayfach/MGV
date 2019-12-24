# Viral detection pipeline

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

<b> Identify circularity (>=20bp direct terminal repeat) </b>  
`python is_circular.py --in_fna $FNA --out_tsv output/circularity.tsv --repeat_length 20`

<b> Estimate completeness </b>  
<b> 1. Download database of complete or circular reference genomes </b>   
Data sources: IMG/VR (XXX circular), UGV (XXX circular), NCBI GenBank (XXX complete)
`wget XXX -O > input/XXX && gunzip input/XXX`  

<b> 2. Use DIAMOND to compute gene sharing and amino acid identity to references </b>
`python compute_aai.py --in_faa $FAA --out_base output/SRS1735492 --threads 10`

<b> 3. Estimate completeness of queries using nearest reference with >80% AAI over >50% genes </b>
