# synteny_plot
Here we use the basics of MCSCAN (python version) to develop a circularized macrosyntenic plot such as the one below:

<img src="https://github.com/user-attachments/assets/0b5040cd-24a2-434b-ae1b-9854ef52672f" alt="Macrosynteny Plot" width="400"/>

Pick 4-6 species of interest. Make sure they have chromosomal level gff files. 

The example we will be using is using 3 cetacean species, 2 terrestrial artiodactyls, and a human to analyze the conserved syntnic regions for the CD300 immune gene family cluster.  

# Find Chromosomes Containing the CD300 Gene Cluster (or your gene of interest)
Phocoena sinus (Vaquita): Chromosome 20 (58,878,645 bp)

Physeter catodon (Sperm Whale): Chromosome 14 (132,988,635 bp)

Homo sapiens (Human): Chromosome 17 (83,257,441 bp)

Camelus dromedarius (Camel): Chromosome 16 (58,371,452 bp)

Bos taurus (Cow): Chromosome 19 (63,449,741 bp)

Balaenoptera musculus (Blue Whale): Chromosome 20 (60,304,989 bp)

# MCSCAN formatting:
- This part is taken directly from https://github.com/tanghaibao/jcvi/wiki/Mcscan-(python-version) --> The scripts are all available there.
  
### Data Preparation:
Extract GFF Data: Convert GFF3 files to BED format to extract mRNA entries.
`python -m jcvi.formats.gff bed --type=mRNA --key=Name Physeter_catodon.ASM283717v2.113.gff3.gz -o sperm_whale_parsed.bed`

### Format FASTA Files: Reformat the CDS FASTA files for easier processing.
`python -m jcvi.formats.fasta format Physeter_catodon.ASM283717v2.cds.all.fa.gz sperm_whale.cds`

### Filter Chromosomes: Filter BED files to only include specific chromosomes.
`awk '$1 == 14' sperm_whale_parsed.bed > sperm_whale_parsed_chrom14.bed`

### Parse CDS: Extract CDS sequences using BED files.
`python parse_cds.py sperm_whale_parsed_chrom14.bed sperm_whale.cds sperm_whale_parsed_chrom14.cds`

### Then run the following command for every species combo: 
`python -m jcvi.compara.catalog ortholog species_a_parsed species_b_parsed --no_strip_names `

(15 combos for CD300 cetacean plot): 

Human - Blue Whale 

Human - Sperm Whale

Human - Vaquita

Human - Camel

Human - Cow

Blue Whale - Sperm Whale

Blue Whale - Vaquita

Blue Whale - Camel

Blue Whale - Cow

Sperm Whale - Vaquita

Sperm Whale - Camel

Sperm Whale - Cow

Vaquita - Camel

Vaquita - Cow 

Camel - Cow 


# Parsing the anchor files: 
- You should now have several anchor files that you need to parse. You need to parse each anchor file from each pairing above to get formatted csv files with synteny and genomic location. Make sure to match to correct bed files.

Run:
  ```Anchors_parsed.ipynb```
For every anchor file. Make sure to match to correct bed files.
  
- You can then grab specific chromosome pairings you want from each parsed csv file: run:``` filter_synteny.ipynb```

# Plot in R
Takes every chromosome specific, parsed csv file and converts it into a circular macrosynteny plot
macrosynteny_plot.r







