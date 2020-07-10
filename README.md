# Consensus Gene Caller
This program attempts to follow the [SEA-PHAGES](https://seaphages.org/) bioinformatic protocol to automate its phage gene annotation. The goal of this program is to reduce the amount of time needed to accurately annotate phage genomes through a customised automation script. The end result produces a spreadsheet with statistics of each gene call to provide information to include or exclude a gene in a final genome annotation. More information about the SEA-PHAGES annotation protocol can be seen [here](https://seaphagesbioinformatics.helpdocsonline.com/home). The algorthim as also designed to be modular, allowing for the additional gene-calling prgrams to be easily included or additional statistical tests.

## Installation
The goal is to make installation of this algorithim as painless as possible, therefore coda installations are preferred. To install conda follow this [link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)


A new environment can be created:
```
conda create -n seaphages -c bioconda -c hcc -c anaconda -c conda-forge prodigal metagene_annotator glimmer fraggenescan hmmer trnascan-se phanotate biopython pandas numpy
conda activate seaphages
git clone https://github.com/ash-bell/Consensus_Gene_Caller.git
```
If you want to run glimmer3 you will need ELPH : Estimated Locations of Pattern Hits from [here](https://cbcb.umd.edu/software/ELPH/). Download and install it as recomended and add it to your PATH. How? [Here](https://opensource.com/article/17/6/set-path-linux)


## Quick start
**Note:** This program can only process one contig at a time. If you have multiple contigs, please run them seperately.
```
./genecaller.py -i phage.fasta -o output_dir -gc all -db pvog pfam -p ~/databases/pvog_database.hmm -f ~/databases/Pfam-A.hmm
./gene_stats.py -i output_dir/gene_clusters.tsv -o outdir2
```

## Programs used
* biopython 1.77
* fraggenescan 1.31
* glimmer 3.02
* hmmer 3.3
* metagene_annotator 1.0
* numpy 1.18.5
* pandas 1.05
* phanotate 2019.08.09
* prodigal 2.6.3
* python 3.8.3
* trnascan-se 2.0.5


## Walkthrough
The program consists of two python scripts, `genecaller.py` and `gene_stats.py`. `genecaller.py` automates the calling of five (I would like to include more in the future) different gencallers prodigal, glimmer, metagene annotator, phanotate and fraggenescan. The output of these are collated and duplicate gene calls removed. The resulting genes can be searched against the Prokaryotic Virus Orthologous Groups (pVOGs) and or Pfam database to find simular genes. It is not required to use all the above listed gencallers. Using `-gc` option, you can specify which genecallers you want. For example `-gc prodigal mga` will just run prodigal and mga. The default is `-gc all`.


By default the search using hidden markov models of the gene call outputs against the pVOG/Pfam database is off, but can be specified with `-db pvog -p <location of pVOG database` and `-db pfam -f <location of pVOG database` or both (see quickstart). This requires the pVOGs and or Pfam database to be setup before hand. The pVOGs database can be downloaded from [here](http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads.html) and the subsequent annotations interpreted from their abbrevations - VOG0002 to "major coat protein" for example. There are plans to include searches against the interpro databases and BLAST against the nr database in the future.

To get the pVOG database working, download the full database (or any subsection), and unzip the folder. Concatinate all the .hmm files and run hmmpress on it to create a database. For example:
```
tar -xvzf AllvogHMMprofiles.tar.gz
cat AllvogHMMprofiles/*.hmm > pVOG
hmmpress pVOG
```

If you want the pFAM-A database, go to ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases and select the lastest release (Pfam33.1/		26/06/2020, 12:02:00 was the lastest for me) and download Pfam-A.hmm.gz. Using `gunzip Pfam-A.hmm.gz` you can unzip the file the database is ready to use.


`gene_stats.py` attempts to closely follow the current annotation protocol popularised by the SEA-PHAGE consortium. You can find more details [here](https://doi.org/10.3390/ijms20143391), which I refer to as the "paper". To use this tool, you just need to specify the output file from `genecaller.py` (the gene_clusters.tsv file) and an output file. The file consists of a spreadsheet with various gene statistics and their scores suggested by the paper to include or exclude a gene from final annotation. Again, this is just a tool to assist in annotation and should not be treated as a final polish product. I will go through each of the methods in the paper and show how I have attempted to impliment a coding alternative.

## Annotation
All the below criteria have been derived from "A Method for Improving the Accuracy and Efficiency of Bacteriophage Genome Annotation" by Salisbury and Tsourkas, 2019 from the International Journal of Molecular Sciences. I would highly recomend giving thier paper a read to understand how I have attemped to automate their methods section.

### 4.1.1. Auto-Annotation Program Program Calls
Annotation is first performed with Glimmer, GeneMark family, Prodigal, PHANOTATE. Glimmer, Prodigal and PHANOTATE is implimented with genecaller.py, the GeneMark family is not. The user could use the webversions of the GeneMark family and integrate their results into the spreadsheet. There are plans to include the automation of the GeneMark family in the future, but for now are ignored. Instead, FragGeneScan and MetaGene Annotator are included. Detection of tRNA genes are not performed in this algorthim with plans to inclde them in the future. (Aragorn, tRNAScan-se ... etc). The algorthim is entirely reliant on the five genecallers to pull out all possible open reading frames (ORFs).

To combined multiple different genecalls together, I have said if a gene has the same start, end and strand(frame), it is the same gene and a representative of that is kept, the rest discarded. Any gene that then has multiple representatives (I call this consenus) is awarded a score of 1 (or True) in the "duplicate" column of the spreadsheet. The paper calls for a point per consensus gene call (for example if prodigal, mga and glimmer all said gene X was in the same place, it would get 3 points). I have only awarded one point because this can unfairly benefit programs that have the same mechanism of calling genes.

The "truncated" column indicates that prodigal or mga thinks the gene is cut in half. For example, it thinks the gene continues off the end of the contig. Should this be the case, you may want to combine the first and last gene together in that order depending on the coding strand (+ or -). If the gene is thought to be truncated it loses a score of 1.

Every genecall provides a score of how "likely" a particular genecaller thinks a gene is correct. Each program calculates its score differently (and annoying not X/100%) so it is diffcult to compare. My alternative is to make each score a fraction of the maximum score achieved. For example if prodgials scores of genes 1-5 are [2, 3, 10, 2, 7] then max(2+3+10+2+7) = 10 (10 is the higest observed score) and the new genescores would be 2/10=0.2, 3/10=0.3, 10/10=1.0 ... etc. This allows a somewhat better comparison of gene scores. (Sorry PHANOTATE). The scores are then add to the final total.

### 4.1.2. Coding Potential
Determining coding potential bioinformatically is diffcult as the output is a graph and subjective to classify. Therefore, some genecallers (prodigal and MGA) provide a ribosomal binding score (RBS). This RBS score is also "normalised" like above and included as the final score. This does make it more likely to include MGA and Prdigal gene calls but RBS is my alternative to coding potential which is difficult to include bioinformatically. rbs_score is also added to the final total.


### 4.1.3. Sequence Similarity Matches
BLASTing or HMMing against a database is timeconsuming and can require setup of large databases to run locally. Therefore, I have made this step optional and default off (although recomended) within this program. The paper suggest doing a pBLAST against the NCBI’s non-redundant (nr) database and include searches of Pfam and Interpro with HMMer. Hits smaller that E-50 are awarded 3 points, between E-50 and E-20 2 points and between E-20 and E-10 1 point. I have implimented the same point system with a HMMer search of the Prokaryotic Virus Orthologous Groups (pVOGs) and Pfam. I plan to include an Interpro HMMer search in the future and maybe a way to integrate a BLAST search against the nr database. (I don't want to do this immediatly due to the nr database size and the time it takes to search, but could provide a script to integrate the data). If hits for the same gene are found against both the Pfam and pVOG database, only the one with the smallest (therefore best) e-value is included in the final scoring system. (Genes with hits from multiple datbases are not double counted)


### 4.1.4. Overlap and Operons
The paper also includes a system to penalise hypothetical genes that overlap with coding genes. It penalises putative genes that overlap by 100 base pair (bp) or more with -4 points, between 100-70 -3 points, 70-40 -2 point and between 40-10 -1 point. I have implimented the same system here but include all genes not just putative genes. This is because calcualing overlap is a timeconsuming process, and better to delete data then need it. I have no included a system to require a 50bp region between forward genes upsteam of reverse genes due to the complexity so do be aware of this. This is because the two genes require space for promoters.


There is also an additional column for operons, and if any overlap is either 1,4 or 8 bps in length it will also append a +1 score to the operon column. This is also included in the final score.


### 4.1.5. Checking ORF Length
Short genes are likely to be coding regions and therefore are penalised (<200bp). Therefore, the paper penalises genes that are between 200 and 150 bps by -1 point, between 150 and 120 -2 points, between 120 and 90 -3 points and smaller than 90bps -4 points. This is implimented within this program and included in the final score.


### 4.1.6. Checking for False Positives
It is up to the user to determine if they wish to keep False positives (genes that are only called for by one program and indicated in the duplicate columm == 0 or False).


### 4.1.7. Decision-Making for Gene Identification
The paper includes a final score which is the sum of all the above criteria. They recomend keeping genes that a score of 3 or greater (don't do this for my program, my scoring system is not identical to theirs). I don't have a recomendation of a score, its up to you to decide.


## Summary
The paper refers to all locations a gene may start as a start codon and uses the criteria of: 1) start codon includes all coding potential (not implimented - "rbs_score" as an alternative); 2) any overlap or operons (implimented as "length_pen" and "operon"); 3) number of prgrams that predict the same start codon (implimented as "duplicate"); 4) high e-value (implimented "e-val_pen"); 5) Shine-Dalgarno score (not implimented - (gene)"score" as an alternative); 6) gene length (implimented as "length_pen").
`dataframe["total_score"] = dataframe[["rbs_score", "score", "duplicate", "length_penality", "operon", "overlap_pen"]].fillna(0).sum(axis=1) - dataframe["truncated"].fillna(0)`

## Future plans
* Include consenus tRNA caller with aragorn, tRNAScan-se
* hmmsearch against interproscan and Pfam
* diamond blast against nr
* Parser for GeneMark family
