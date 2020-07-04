#!/usr/bin/env python3

import subprocess
import python_colours as clr
import argparse
import pandas as pd
import logging
import os
import numpy as np

logging.basicConfig(level=logging.INFO, filename='logfile', filemode='w', format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger()


def my_args():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Consenus Genecalling of a FASTA file")

    # Positional mandatory arguments
    parser.add_argument("--infile", "-i", dest="infile", required=True, help="FASTA file in nuclotide format")
    parser.add_argument("--outdir", "-o", dest="output_folder", required=True, help="Output folder where the files will be created or generated if it doesn't exist")
    parser.add_argument("--genecallers", "-gc", dest="genecallers",  nargs="+", default=["all"], choices=["prodigal", "mga", "glimmer", "fraggenescan", "phanotate", "all"],
                        help="Select genecaller(s)")
    # Parse arguments
    args = parser.parse_args()

    return args


def create_output_dir(directory_name):
    try:
        os.mkdir(directory_name)
        logger.debug('Created output folder at {}'.format(args.output_folder))
    except OSError:
        logger.debug('Output folder at {} already exists'.format(args.output_folder))


def runProdigal(infile):
    '''
    Run the prodigal gene caller on the input FASTA file and output the gene calls in a table.
    '''
    print(f"{clr.colours.BGreen}Calling genes with prodigal{clr.colours.Colour_Off}")
    cmd = f"""prodigal -i {infile} -f gff -p meta -q -m -g 11 -o {args.output_folder}/prodigal.txt"""
    stdout, stderr = execute(cmd)
    prodigal = pd.read_csv(f"{args.output_folder}/prodigal.txt", sep="\t|;", comment='#', engine='python', header=None,
                            usecols=[0,3,4,5,6,8,9,18], names=["contig", "start", "end", "score", "strand", "gene", "truncated", "rbs_score"])
    prodigal["rbs_score"] = prodigal["rbs_score"].str.split('=').str[1]
    prodigal["rbs_score"] = pd.to_numeric(prodigal["rbs_score"])
    prodigal["rbs_score"] = prodigal["rbs_score"] / max(prodigal["rbs_score"])
    prodigal["truncated"] = prodigal["truncated"].str.split('=').str[1]
    prodigal["truncated"] = pd.to_numeric(prodigal["truncated"])
    prodigal.loc[prodigal["truncated"] > 0, "truncated"] = "True"
    prodigal.loc[prodigal["truncated"] == 00, "truncated"] = "False"
    prodigal["gene"] = "prodigal_"+prodigal["gene"].str.split('=').str[1]
    prodigal["score"] = prodigal["score"] / max(prodigal["score"])
    prodigal.to_csv(f"{args.output_folder}/prodigal.out", index=False)


def runMGA(infile):
    '''
    Run runMGA
    '''
    print(f"{clr.colours.BGreen}Calling genes with MetaGene Annotator{clr.colours.Colour_Off}")
    cmd = (f"""mga {infile} -m > {args.output_folder}/mga.txt""")
    stdout, stderr = execute(cmd)
    mga = pd.read_csv(f"{args.output_folder}/mga.txt", sep="\t", comment='#', engine='python', header=None,
                        usecols=[0,1,2,3,5,6,10], names=["gene", "start", "end", "strand", "truncated", "score", "rbs_score"])
    mga["gene"] = "MGA_"+mga["gene"]
    with open(f"{args.output_folder}/mga.txt", "r") as handle:
        first_line = handle.readline()
    mga["contig"] = first_line.split()[1]
    mga["truncated"] = pd.to_numeric(mga["truncated"])
    mga.loc[mga["truncated"] <= 10, "truncated"] = "True"
    mga.loc[mga["truncated"] == 11, "truncated"] = "False"
    mga["score"] = mga["score"] / max(mga["score"])
    mga["rbs_score"].replace("-", np.nan, inplace=True)
    mga["rbs_score"] = pd.to_numeric(mga["rbs_score"])
    mga["rbs_score"] = mga["rbs_score"] / max(mga["rbs_score"])

    mga.to_csv(f"{args.output_folder}/mga.out", index=False)


def runGlimmer(infile):
    '''
    '''
    print(f"{clr.colours.BGreen}Calling genes with Glimmer3{clr.colours.Colour_Off}")
    cmd = f"""
            mkdir {args.output_folder}/GLIMMER;
            long-orfs -n -t 1.15 {infile} {args.output_folder}/GLIMMER/run3.longorfs;
            extract -t {infile} {args.output_folder}/GLIMMER/run3.longorfs > {args.output_folder}/GLIMMER/run3.train;
            build-icm -r {args.output_folder}/GLIMMER/run3.icm < {args.output_folder}/GLIMMER/run3.train;
            glimmer3 -o50 -g110 -t30 {infile} {args.output_folder}/GLIMMER/run3.icm {args.output_folder}/GLIMMER/run3.run1;
            tail -n +2 {args.output_folder}/GLIMMER/run3.run1.predict > {args.output_folder}/GLIMMER/run3.coords;
            upstream-coords.awk 25 0 {args.output_folder}/GLIMMER/run3.coords | extract {infile} - > {args.output_folder}/GLIMMER/run3.upstream;
            elph {args.output_folder}/GLIMMER/run3.upstream LEN=6 | get-motif-counts.awk > {args.output_folder}/GLIMMER/run3.motif;
            startuse="$(start-codon-distrib -3 {infile} {args.output_folder}/GLIMMER/run3.coords)";
            glimmer3 -o50 -g110 -t30 -b {args.output_folder}/GLIMMER/run3.motif -P $startuse {infile} {args.output_folder}/GLIMMER/run3.icm {args.output_folder}/glimmer;
            rm -r {args.output_folder}/GLIMMER {args.output_folder}/glimmer.detail
            """
    stdout, stderr = execute(cmd)
    glimmer = pd.read_csv(f"{args.output_folder}/glimmer.predict", sep="\s+", comment=">", engine="python", header=None, names=["gene", "start", "end", "strand", "score"])
    glimmer["gene"] = "glimmer_"+glimmer["gene"]
    with open(f"{args.output_folder}/glimmer.predict", "r") as handle:
        first_line = handle.readline()
    glimmer["contig"] = first_line.split()[0].split(">")[1]
    glimmer.loc[glimmer["strand"] < 0, "strand"] = "-"
    glimmer.loc[glimmer["strand"] > 0, "strand"] = "+"
    glimmer["score"] = glimmer["score"] / max(glimmer["score"])
    glimmer["rbs_score"] = np.nan
    glimmer["truncated"] = np.nan
    glimmer.to_csv(f"{args.output_folder}/glimmer.out", index=False)


def runFragGeneScan(infile):
    '''
    '''
    print(f"{clr.colours.BGreen}Calling genes with FragGeneScan{clr.colours.Colour_Off}")
    cmd = f"""run_FragGeneScan.pl -genome={infile} -out={args.output_folder}/FGS -complete=1 -train=complete"""
    stdout, stderr = execute(cmd)
    data = pd.read_csv(f"{args.output_folder}/FGS.gff", sep="\t|;", comment="#", engine="python", header=None,
                        usecols=[0,3,4,6,8], names=["contig", "start", "end", "strand", "gene"])
    df = pd.read_csv(f"{args.output_folder}/FGS.out", sep="\t", comment=">", engine="python", header=None,
                        usecols=[0,1,2,4], names=["start", "end", "strand", "score"])
    FGS = pd.merge(data, df, on=["start", "end", "strand"], how="outer")
    FGS["rbs_score"] = np.nan
    FGS["truncated"] = np.nan
    FGS["gene"] = "FGS_"+FGS["gene"]
    FGS["score"] = FGS["score"] / max(FGS["score"])
    FGS.to_csv(f"{args.output_folder}/fgs.out", index=False)


def runPHANOTATE(infile):
    '''
    '''
    print(f"{clr.colours.BGreen}Calling genes with PHANOTATE{clr.colours.Colour_Off}")
    cmd = f"""phanotate.py -o {args.output_folder}/PHANOTATE.out -f tabular {infile}"""
    stdout, stderr = execute(cmd)
    PHANOTATE = pd.read_csv(f"{args.output_folder}/PHANOTATE.out", sep="\t", comment="#", engine="python",
                            header=None, names=["start", "end", "strand", "contig", "score"], index_col=False)
    PHANOTATE["rbs_score"] = np.nan
    PHANOTATE["truncated"] = np.nan
    PHANOTATE["gene"] = "PHANOTATE_"+PHANOTATE.index.astype(str)
    PHANOTATE["score"] = PHANOTATE["score"] / min(PHANOTATE["score"])
    PHANOTATE.to_csv(f"{args.output_folder}/phan.out", index=False)


def contig_length(infile):
    '''
    get the contig length incase gene is truncated
    '''
    cmd = f"""grep -v ">" {infile} | wc -c > {args.output_folder}/contig_length.txt"""
    stdout, stderr = execute(cmd)


def combineGeneClusters():
    '''
    '''
    files = ["prodigal.out", "mga.out", "glimmer.out", "fgs.out", "phan.out"]
    list = []
    for f in files:
        if os.path.isfile(f"{args.output_folder}/{f}"):
            df = pd.read_csv(f"{args.output_folder}/{f}")
            list.append(df)
            os.remove(f"{args.output_folder}/{f}")
        else:
            pass

    frame = pd.concat(list, sort=True, ignore_index=True)
    frame["duplicate"] = frame.duplicated(subset=["contig", "start", "end", "strand"], keep=False)
    frame = frame.drop_duplicates(["contig", "start", "end", "strand"]).sort_values(by=["start"]).reset_index(drop=True)

    with open(f"{args.output_folder}/contig_length.txt", "r") as handle:
        contig_length = handle.readline()
    frame["contig_length"] = contig_length

    frame.to_csv(f"{args.output_folder}/gene_clusters.tsv", index=None)


def execute(bash):
    logger.info('Executing {}'.format(bash))
    process = subprocess.Popen(bash, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (stdout, stderr) = process.communicate()

    logger.debug(stdout.decode("utf-8"), stderr.decode("utf-8"))
    return stdout, stderr


def main(args):
    create_output_dir(args.output_folder)
    contig_length(args.infile)

    if "all" in args.genecallers:
        runProdigal(args.infile,)
        runMGA(args.infile,)
        runGlimmer(args.infile,)
        runFragGeneScan(args.infile,)
        runPHANOTATE(args.infile,)
    if "prodigal" in args.genecallers:
        runProdigal(args.infile,)
    if "mga" in args.genecallers:
        runMGA(args.infile,)
    if "glimmer" in args.genecallers:
        runGlimmer(args.infile,)
    if "fraggenescan" in args.genecallers:
        runFragGeneScan(args.infile,)
    if "phanotate" in args.genecallers:
        runPHANOTATE(args.infile,)

    combineGeneClusters()


if __name__ == "__main__":
    args = my_args()
    main(args)
