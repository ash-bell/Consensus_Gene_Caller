#!/usr/bin/env python3

import subprocess
import python_colours as clr
import argparse
import pandas as pd
import logging
import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
    parser.add_argument("--trna", dest="trna", nargs="+", default=["all"], choices=["all", "aragorn", "trnascan"], help="Select tRNA Scanner")
    parser.add_argument("--database", "-db", dest="database", default=["None"], nargs="+", choices=["pvog", "pfam"], help="Select Database")
    parser.add_argument("--pvog", "-p", dest="pvog", default=["None"], help="pVOG HMM database file")
    parser.add_argument("--pfam", "-f", dest="pfam", default=["None"], help="Pfam HMM database file")
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
    if os.path.isfile(f"{args.output_folder}/prodigal.out"):
        print(f"{clr.colours.BWhite}Prodigal output file {args.output_folder}/prodigal.out exist, skipping prodigal{clr.colours.Colour_Off}")
        logger.info(f'Prodigal output file {args.output_folder}/prodigal.out exist, skipping prodigal')
    else:
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
        os.remove(f"{args.output_folder}/prodigal.txt")


def runMGA(infile):
    '''
    Run runMGA
    '''
    if os.path.isfile(f"{args.output_folder}/mga.out"):
        print(f"{clr.colours.BWhite}MetaGene Annotator output file {args.output_folder}/mga.out exist, skipping MGA{clr.colours.Colour_Off}")
        logger.info(f'MGA output file {args.output_folder}/mga.out exist, skipping MGA')
    else:
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
        os.remove(f"{args.output_folder}/mga.txt")


def runGlimmer(infile):
    '''
    '''
    if os.path.isfile(f"{args.output_folder}/glimmer.out"):
        print(f"{clr.colours.BWhite}Glimmer output file {args.output_folder}/glimmer.out exist, skipping glimmer3{clr.colours.Colour_Off}")
        logger.info(f'Glimmer output file {args.output_folder}/glimmer.out exist, skipping glimmer3')
    else:
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
        os.remove(f"{args.output_folder}/glimmer.predict")


def runFragGeneScan(infile):
    '''
    '''
    if os.path.isfile(f"{args.output_folder}/fgs.out"):
        print(f"{clr.colours.BWhite}FragGeneScane output file {args.output_folder}/fgs.out exist, skipping FGS{clr.colours.Colour_Off}")
        logger.info(f'FragGeneScane output file {args.output_folder}/fgs.out exist, skipping FGS')
    else:
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
        os.remove(f"{args.output_folder}/FGS.faa")
        os.remove(f"{args.output_folder}/FGS.ffn")
        os.remove(f"{args.output_folder}/FGS.gff")
        os.remove(f"{args.output_folder}/FGS.out")


def runPHANOTATE(infile):
    '''
    '''
    if os.path.isfile(f"{args.output_folder}/phan.out"):
        print(f"{clr.colours.BWhite}PHANOTATE output file {args.output_folder}/phan.out exist, skipping PHANOTATE{clr.colours.Colour_Off}")
        logger.info(f'PHANOTATE output file {args.output_folder}/phan.out exist, skipping PHANOTATE')
    else:
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
        os.remove(f"{args.output_folder}/PHANOTATE.out")


def contig_length(infile):
    '''
    get the contig length incase gene is truncated
    '''
    cmd = f"""grep -v ">" {infile} | wc -c > {args.output_folder}/contig_length.txt"""
    stdout, stderr = execute(cmd)


def combineGeneClusters():
    '''
    '''
    print(f"{clr.colours.BGreen}Combining simular gene clusters based on start/end location{clr.colours.Colour_Off}")
    files = ["prodigal.out", "mga.out", "glimmer.out", "fgs.out", "phan.out"]
    list = []
    for f in files:
        if os.path.isfile(f"{args.output_folder}/{f}"):
            df = pd.read_csv(f"{args.output_folder}/{f}")
            list.append(df)
        else:
            pass

    df = pd.concat(list, sort=True, ignore_index=True)
    df["duplicate"] = df.duplicated(subset=["contig", "start", "end", "strand"], keep=False)
    df = df.drop_duplicates(["contig", "start", "end", "strand"]).sort_values(by=["start"]).reset_index(drop=True)

    with open(f"{args.output_folder}/contig_length.txt", "r") as handle:
        contig_length = handle.readline()
    df["contig_length"] = int(contig_length)
    df.to_csv(f"{args.output_folder}/gene_clusters.tsv", sep="\t", index=None)
    os.remove(f"{args.output_folder}/contig_length.txt")
    return(df)


def runAragorn(infile):
    if os.path.isfile(f"{args.output_folder}/aragorn.out"):
        print(f"{clr.colours.BWhite}aragorn output file {args.output_folder}/aragorn.out exist, skipping aragorn{clr.colours.Colour_Off}")
    else:
        print(f"{clr.colours.BGreen}Running ARAGORN a tRNA finder{clr.colours.Colour_Off}")
        cmd = f"""aragorn -gc11 -w -Q {infile} | awk '{{print $2, $3}}' > {args.output_folder}/aragorn.out"""
        stdout, stderr = execute(cmd)


def runtRNAscan(infile):
    if os.path.isfile(f"{args.output_folder}/tRNAscan.out"):
        print(f"{clr.colours.BWhite}tRNAscan output file {args.output_folder}/tRNAscan.out exist, skipping tRNAscan{clr.colours.Colour_Off}")
    else:
        print(f"{clr.colours.BGreen}Running tRNAscan-SE 2.0 a tRNA finder{clr.colours.Colour_Off}")
        cmd = f"tRNAscan-SE -B -o {args.output_folder}/tRNAscan.out -Q {infile}"
        stdout, stderr = execute(cmd)


def combinetRNAcalls():
    print(f"{clr.colours.BGreen}Combining tRNA calls{clr.colours.Colour_Off}")
    aragorn = pd.read_csv(f"{args.output_folder}/aragorn.out", sep="\s+", skiprows=2, header=None, names=["aragorn", "coding"])
    aragorn[["start", "end"]] = aragorn["coding"].str.split(",", expand=True).replace(regex=[r'c\[', '\[', '\]'], value='').astype('int64')
    aragorn.drop(["coding"], axis=1, inplace=True)

    tRNAscan = pd.read_csv(f"{args.output_folder}/tRNAscan.out", sep="\t", skiprows=3, header=None, usecols=[2, 3, 4], names=["start", "end", "trnascan"])

    df = pd.merge(aragorn, tRNAscan, on=["start", "end"], how="outer", indicator=True)
    new_df = df[df["_merge"] != "right_only"]
    data = df[df["_merge"] == "right_only"]
    flipped = data.rename(columns={"start": "end", "end": "start"})
    consensus = pd.merge(new_df, flipped, on=["start", "end"], how="outer")
    consensus.drop(["_merge_x", "_merge_y"], axis=1, inplace=True)
    consensus.to_csv(f"{args.output_folder}/tRNA_calls.tsv", sep="\t")
    print(f"""{clr.colours.BBlue}See file {args.output_folder}/tRNA_calls.tsv for details{clr.colours.Colour_Off}""")

def getSeq(df):
    print(f"{clr.colours.BGreen}Getting consenus gene calls{clr.colours.Colour_Off}")
    record = SeqIO.read(f"{args.infile}", "fasta")
    sequences = []
    for i in range(len(df)):
        if df["strand"][i] == "+":
            if df["start"][i] > df["end"][i]:
                sequences.append(SeqRecord(record.seq[df["start"][i]-1:len(record)]+record.seq[0:df["end"][i]], id=df["gene"][i], description=""))
            else:
                sequences.append(SeqRecord(record.seq[df["start"][i]-1:df["end"][i]], id=df["gene"][i], description=""))
        if df["strand"][i] == "-":
            if df["end"][i] > df["start"][i]:
                sequences.append(SeqRecord((record.seq[df["end"][i]-1:len(record)]+record.seq[0:df["start"][i]]).reverse_complement(), id=df["gene"][i], description=""))
            else:
                sequences.append(SeqRecord(record.seq[df["end"][i]-1:df["start"][i]].reverse_complement(), id=df["gene"][i], description=""))
    with open(f"{args.output_folder}/all_gene_calls.fna", "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
    return(sequences)


def translateSeq(seq_record):
    print(f"{clr.colours.BGreen}Translating coding regions{clr.colours.Colour_Off}")
    for record in seq_record:
        remainder = len(record.seq) % 3
        if remainder == 0:
            record.seq = record.seq.translate(table=11)
        else:
            record.seq = (record.seq + Seq('N' * (3 - remainder))).translate(table=11)
    with open(f"{args.output_folder}/all_gene_calls.faa", "w") as output_handle:
        SeqIO.write(seq_record, output_handle, "fasta")


def runpVOGsearch():
    if args.pvog == ["None"]:
        print(f"{clr.colours.BRed}no pVOG database provided, HMM search against pVOG database will not be run{clr.colours.Colour_Off}")
    elif os.path.isfile(f"{args.output_folder}/pVOG_hits.tbl"):
        pVOGS = pd.read_csv(f"{args.output_folder}/pVOG_hits.tbl", sep="\s+", comment="#", engine='python', header=None,
                            usecols=[0, 2, 4], names=["gene", "pVOG_hits", "pVOG_e-val"])
        pVOGS.sort_values("pVOG_e-val", inplace=True)
        pVOGS.drop_duplicates(subset="gene", keep="first", inplace=True)
        df = pd.read_csv(f"{args.output_folder}/gene_clusters.tsv", sep="\t")
        data = df.merge(pVOGS, on="gene", how="outer")
        data.to_csv(f"{args.output_folder}/gene_clusters.tsv", sep="\t")
    elif os.path.isfile(f"{args.pvog}"):
        print(f"{clr.colours.BGreen}Searching protein coding regions against the pVOG database {clr.colours.Colour_Off}")
        cmd = f"""hmmsearch --tblout {args.output_folder}/pVOG_hits.tbl --notextw -E 10 {args.pvog} {args.output_folder}/all_gene_calls.faa"""
        stdout, stderr = execute(cmd)
        pVOGS = pd.read_csv(f"{args.output_folder}/pVOG_hits.tbl", sep="\s+", comment="#", engine='python', header=None,
                            usecols=[0, 2, 4], names=["gene", "pVOG_hits", "pVOG_e-val"])
        pVOGS.sort_values("pVOG_e-val", inplace=True)
        pVOGS.drop_duplicates(subset="gene", keep="first", inplace=True)
        df = pd.read_csv(f"{args.output_folder}/gene_clusters.tsv", sep="\t")
        data = df.merge(pVOGS, on="gene", how="outer")
        data.to_csv(f"{args.output_folder}/gene_clusters.tsv", sep="\t")
    else:
        print(f"""{clr.colours.BRed}'{args.pvog}' should be a file path and is not. Please specify the pVOG database filepath{clr.colours.Colour_Off}""")


def runPfamsearch():
    if args.pfam == ["None"]:
        print(f"{clr.colours.BRed}no Pfam database provided, HMM search against Pfam database will not be run{clr.colours.Colour_Off}")
    elif os.path.isfile(f"{args.output_folder}/Pfam_hits.tbl"):
        Pfam = pd.read_csv(f"{args.output_folder}/Pfam_hits.tbl", sep="\s+", comment="#", engine='python', header=None,
                            usecols=[0, 2, 4], names=["gene", "Pfam_hit", "Pfam_e-val"])
        Pfam.sort_values("Pfam_e-val", inplace=True)
        Pfam.drop_duplicates(subset="gene", keep="first", inplace=True)
        df = pd.read_csv(f"{args.output_folder}/gene_clusters.tsv", sep="\t")
        data = df.merge(Pfam, on="gene", how="outer")
        data.to_csv(f"{args.output_folder}/gene_clusters.tsv", sep="\t")
    elif os.path.isfile(f"{args.pfam}"):
        print(f"{clr.colours.BGreen}Searching protein coding regions against the Pfam-A database {clr.colours.Colour_Off}")
        cmd = f"""hmmsearch --tblout {args.output_folder}/Pfam_hits.tbl --notextw -E 10 {args.pfam} {args.output_folder}/all_gene_calls.faa"""
        stdout, stderr = execute(cmd)
        Pfam = pd.read_csv(f"{args.output_folder}/Pfam_hits.tbl", sep="\s+", comment="#", engine='python', header=None,
                            usecols=[0, 2, 4], names=["gene", "Pfam_hit", "Pfam_e-val"])
        Pfam.sort_values("Pfam_e-val", inplace=True)
        Pfam.drop_duplicates(subset="gene", keep="first", inplace=True)
        df = pd.read_csv(f"{args.output_folder}/gene_clusters.tsv", sep="\t")
        data = df.merge(Pfam, on="gene", how="outer")
        data.to_csv(f"{args.output_folder}/gene_clusters.tsv", sep="\t")
    else:
        print(f"""{clr.colours.BRed}'{args.pfam}' should be a file path and is not. Please specify the Pfam database filepath{clr.colours.Colour_Off}""")


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

    if "all" in args.trna:
        runAragorn(args.infile)
        runtRNAscan(args.infile)
        if os.path.getsize(f"{args.output_folder}/tRNAscan.out") == 0:
            print(f"{clr.colours.BWhite}tRNAscan-SE 2.0 did not find any tRNAs{clr.colours.Colour_Off}")
        if os.path.getsize(f"{args.output_folder}/aragorn.out") == 14:
            print(f"{clr.colours.BWhite}aragorn did not find any tRNAs{clr.colours.Colour_Off}")
        else:
            combinetRNAcalls()
    if "aragorn" in args.trna:
        runAragorn(args.infile)
        print(f"""{clr.colours.BBlue}See file {args.output_folder}/aragorn.out for details{clr.colours.Colour_Off}""")
    if "trnascan" in args.trna:
        runtRNAscan(args.infile)
        print(f"""{clr.colours.BBlue}See file {args.output_folder}/tRNAscan.out for details{clr.colours.Colour_Off}""")

    if args.database == ["None"]:
        print(f"{clr.colours.BWhite}No databases indicated to search{clr.colours.Colour_Off}")

    if "pvog" in args.database:
        if os.path.isfile(f"{args.output_folder}/all_gene_calls.faa"):
            runpVOGsearch()
        else:
            translateSeq(getSeq(combineGeneClusters()))
            runpVOGsearch()
    if "pfam" in args.database:
        if os.path.isfile(f"{args.output_folder}/all_gene_calls.faa"):
            runPfamsearch()
        else:
            translateSeq(getSeq(combineGeneClusters()))
            runPfamsearch()

    print(f"{clr.colours.BBlue}See file: {args.output_folder}/gene_clusters.tsv {clr.colours.Colour_Off}")


if __name__ == "__main__":
    args = my_args()
    main(args)
