# synteny_plot
Here we use the basics of MCSCAN (python version) to develop a circularized macrosyntenic plot such as the one below:

<img src="https://github.com/user-attachments/assets/0b5040cd-24a2-434b-ae1b-9854ef52672f" alt="Macrosynteny Plot" width="400"/>

Pick 4-6 species of interest. Make sure they have chromosomal level gff files. 

The example we will be using is analyzing the conserved syntney for cetacean CD300 gene family. 

# Find Chromosomes Containing the CD300 Gene Cluster (or your gene of interest)
Phocoena sinus (Vaquita): Chromosome 20 (58,878,645 bp)
Physeter catodon (Sperm Whale): Chromosome 14 (132,988,635 bp)
Homo sapiens (Human): Chromosome 17 (83,257,441 bp)
Camelus dromedarius (Camel): Chromosome 16 (58,371,452 bp)
Bos taurus (Cow): Chromosome 19 (63,449,741 bp)
Balaenoptera musculus (Blue Whale): Chromosome 20 (60,304,989 bp)
Analysis Steps

Data Preparation:
Extract GFF Data: Convert GFF3 files to BED format to extract mRNA entries.

`python -m jcvi.formats.gff bed --type=mRNA --key=Name Physeter_catodon.ASM283717v2.113.gff3.gz -o sperm_whale_parsed.bed`

Format FASTA Files: Reformat the CDS FASTA files for easier processing.

python -m jcvi.formats.fasta format Physeter_catodon.ASM283717v2.cds.all.fa.gz sperm_whale.cds
Filter Chromosomes (Optional): If necessary, filter BED files to only include specific chromosomes.

awk '$1 == 14' sperm_whale_parsed.bed > sperm_whale_parsed_chrom14.bed
Parse CDS: Extract CDS sequences using BED files.

python parse_cds.py sperm_whale_parsed.bed sperm_whale.cds sperm_whale_parsed.cds
For specific chromosomes:


python parse_cds.py sperm_whale_parsed_chrom14.bed sperm_whale.cds sperm_whale_parsed_chrom14.cds
