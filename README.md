# Consensus Gene Caller
This program attempts to follow the SEA-PHAGES bioinformatic protocol to automate its phage gene annotation. The goal of this program is to reduce the amount of time needed to accurately annotate phage genomes through a customised automation script. The end result produces a spreadsheet with statistics of each gene call to provide information to include or exclude a gene in a final genome annotation. More information about the SEA-PHAGES annotation protocol can be seen (here)[https://seaphagesbioinformatics.helpdocsonline.com/home] or reading the PDF included. The algorthim as also designed to be modular, allowing for the additional gene-calling prgrams to be easily included or additional statistical tests.

## Installation
The goal is to make installation of this algorithim as painless as possible, therefore coda installations are preferred. To install conda follow this (link)[https://docs.conda.io/projects/conda/en/latest/user-guide/install/]


A new environment can be created:
conda create -n seaphages -c bioconda -c hcc prodigal metagene_annotator glimmer fraggenescan hmmer trnascan-se phanotate

install elph


## Quick start
```
./genecaller.py -i phage.fasta -o output_dir -gc all -db pvog -p ~/databases/pvog_database.hmm
./gene_stats.py -i output_dir/gene_clusters.tsv -o outdir2
```

## Programs used
biopython 1.77
fraggenescan 1.31
glimmer 3.02
hmmer 3.3
metagene_annotator 1.0
numpy 1.18.5
pandas 1.05
phanotate 2019.08.09
prodigal 2.6.3
python 3.8.3
trnascan-se 2.0.5


## Walkthrough
The program consists of two python scripts, `genecaller.py` and `gene_stats.py`. `genecaller.py` automates the calling of five (I would like to include more in the future) different gencallers prodigal, glimmer, metagene annotator, phanotate and fraggenescan. The output of these are collated and duplicate gene calls removed. The resulting genes can be searched against the Prokaryotic Virus Orthologous Groups (pVOGs) database to find simular genes. It is not required to use all the above listed gencallers. Using `-gc` option, you can specify which genecallers you want. For example `-gc prodigal mga` will just run prodigal and mga. The default is `-gc all`.


By default the search using hidden markov models of the gene call outputs against the pVOG database is off, but can be specified with `-db pvog -p <location of pVOG database`. This requires the pVOGs database to be setup before hand. The pVOGs database can be downloaded from (here)[http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads.html] and the subsequent annotations interpreted from their abbrevations - VOG0002 to "major coat protein" for example. There are plans to include searches against the interpro database and BLAST against the nr database in the future.

`gene_stats.py` attempts to closely follow the current annotation protocol popularised by the SEA-PHAGE consortium. You can find more details (here)[https://doi.org/10.3390/ijms20143391]. 





## Future plans
Include consenus tRNA caller with aragorn, tRNAScan-se
hmmsearch against interproscan
diamond blast against nr
