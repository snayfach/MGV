# Viral detection pipeline

<b>  Install dependencies </b>  
Prodigal v2.6.3   
HMMER v3.1b2  
VirFinder v1.1  

You should be able to run the following from the command line:  
`prodigal`  
`hmmsearch`

And import 'VirFinder' from the R command line:  
`library(VirFinder)`  

<b> Move to subdirectory in repository for viral detection  </b>  
`cd /path/to/UGV/0_viral_detect`  

<b> Download & decompress reference databases  </b>  
`wget -O input/imgvr.hmm.gz https://img.jgi.doe.gov/virus/doc/final_list.hmms.gz`  
`wget -O input/pfam.hmm.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz`
`gunzip input/imgvr.hmm.gz`  
`gunzip input/pfam.hmm.gz`  

<b> Locate input metagenomic assembly </b>  
An example file has been provided for sample_id 'SRS1735492' from project 'JGI' which was analyzed in the manuscript:  
FNA=input/SRS1735492.fna  

<b> Call viral genes </b>  
FAA=input/SRS1735492.faa  
FFN=input/SRS1735492.ffn  
GFF=input/SRS1735492.gff  
`prodigal -i $FNA -a $FAA -d $FFN -p meta > $GFF`

<b> Run HMMER on imgvr and pfam databases </b>  
`hmmsearch -Z 1 --cpu 32 --noali --tblout output/imgvr.out input/imgvr.hmm $FAA`
`hmmsearch -Z 1 --cut_tc --cpu 32 --noali  --tblout output/pfam.out input/pfam.hmm $FAA`

<b> Count genes hitting viral and microbial marker genes </b>  
Each gene assigned according to its best hit E-value <1e-10. Exclude HMMs listed in exclude_hmms directory.  

`python count_hmm_hits.py $FNA $FAA output/imgvr.out output/pfam.out > output/hmm_hits.tsv`

<b> Run VirFinder </b>  
`Rscript virfinder.R $FNA output/virfinder.tsv`  

<b> Quantify strand switch rate of genes </b>  
`python strand_switch.py $FNA $FAA > output/strand_switch.tsv`  

<b> Create master table of sequence features </b>  
`python master_table.py output/hmm_hits.tsv output/virfinder.tsv output/strand_switch.tsv > output/master_table.tsv`  

<b> Predict viral contigs </b>  
```
python viral_classify.py \
--in_features output/master_table.tsv \
--out_features output/SRS1735492.tsv \
--in_seqs_base input/SRS1735492 \
--out_seqs_base output/SRS1735492
```





