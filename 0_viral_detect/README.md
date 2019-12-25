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

<b> Download & decompress reference databases of HMMs  </b>  
`wget -O input/imgvr.hmm.gz https://img.jgi.doe.gov/virus/doc/final_list.hmms.gz`  
`wget -O input/pfam.hmm.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz`
`gunzip input/imgvr.hmm.gz`  
`gunzip input/pfam.hmm.gz`  

<b> Locate input metagenomic assembly </b>  
An example file 'input/SRS1735492.fna' has been provided. This represents metagenomic contigs for sample_id 'SRS1735492' from project 'JGI'  

<b> Call viral genes </b>  
`prodigal -i input/SRS1735492.fna -a input/SRS1735492.faa -p meta`

<b> Run HMMER on imgvr and pfam databases </b>  
`hmmsearch -Z 1 --cpu 32 --noali --tblout output/imgvr.out input/imgvr.hmm input/SRS1735492.faa`
`hmmsearch -Z 1 --cut_tc --cpu 32 --noali  --tblout output/pfam.out input/pfam.hmm input/SRS1735492.faa`

Here the `-Z 1` flag is specified to make E-values comparable between databases and between samples. The IMG/VR database contains viral genes while the Pfam database contains non-viral genes.

<b> Count genes hitting viral and microbial marker genes </b>  
 `python count_hmm_hits.py input/SRS1735492.fna input/SRS1735492.faa output/imgvr.out output/pfam.out > output/hmm_hits.tsv`
 
Each gene assigned according to its best hit with E-value <1e-10. We found that several of the HMMs from IMG/VR HMMs are commonly found in non-viral genomes and several of the HMMs from Pfam are commonly found in viral genomes. Those HMMs are excluded from the analysis.

<b> Run VirFinder </b>  
`Rscript virfinder.R $FNA output/virfinder.tsv`  

VirFinder scores each contig using a machine-learning algorithm based on kmers  
https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0283-5

<b> Quantify strand switch rate of genes </b>  
`python strand_switch.py input/SRS1735492.fna input/SRS1735492.faa > output/strand_switch.tsv`  

Here the code scans the proteins from each genome in genomic order. It counts the nuber of strand switches (+ to - or - to +) and divides by the total number of genes.

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

The code now classifies each sequence as viral based on its combination of viral signatures. These decisions are based on a combination of rules that are specific to different fragment lengths. Each rule is a combination of 'ANDs' (ex: a contig must have >20% viral genes and <5% of non-viral genes) and rules are combined with an 'OR' operator (ex: a contig is viral if it satisfies rule 1 OR rule 2 or rule 3). Up to 5 rule combinations are used per contig length. These rules are listed in the table 'classification_rules.tsv'




