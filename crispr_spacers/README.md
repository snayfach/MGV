# CRISPR-spacer identification

<b> Identify CRISPR arrays using CRT and PILERCR </b>  
python identify_crispr.py -i example/GUT_GENOME147678.fna -o out -p crt

The programs CRT and PILERCR are located in the `bin` directory
The program first splits your input sequence into chunks of 50 Mbp. This is done to limit RAM usage. Not an an issue for small genomes, but can be one for large metagenomes.
Each sequence is padded with 50 Ns, which helps to find CRISPR arrays at contig boundaries
CRT v1.2 and PILER v1.06 are then run on each chunk of data
The results are parsed and written to the output directory

<b>Merge overlapping CRISPR arrays identified using CRT and PILERCR </b>  
python merge_crispr.py out/crt out/pilercr out/merged

