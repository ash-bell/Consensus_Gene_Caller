#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import logging
import numpy as np
import python_colours as clr


logging.basicConfig(level=logging.INFO, filename='logfile_2', filemode='w', format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger()


def my_args():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Calculate Gene Statistics")

    # Positional mandatory arguments
    parser.add_argument("--infile", "-i", dest="infile", required=True, help="File generated from genecaller.py (gene_clusters.tsv)")
    parser.add_argument("--outdir", "-o", dest="output_folder", required=True, help="Output folder where the files will be created or generated if it doesn't exist")

    # Parse arguments
    args = parser.parse_args()

    return args


def create_output_dir(directory_name):
    try:
        os.mkdir(directory_name)
        logger.debug('Created output folder at {}'.format(args.output_folder))
    except OSError:
        logger.debug('Output folder at {} already exists'.format(args.output_folder))


def gene_length(row):
    if row["strand"] == "+":
        if row["start"] > row["end"]:
            return(row["contig_length"] - row["start"] + row["end"])
        else:
            return(row["end"] - row["start"])
    if row["strand"] == "-":
        if row["end"] > row["start"]:
            return(row["contig_length"] - row["end"] + row["start"])
        else:
            return(row["start"] - row["end"])
    return(np.nan)


def length_pen(dataframe):
    '''
    1) Calculate gene length. This allows for standardisation of how Gene Length is defined as Glimmer's is different (3bps shorter then GeneMark Family for instance).
    2) Create a penality system for short genes.
        Scoring
        -------
        If Gene length is > 200 bps, 0 points
        If Gene length is between 200 and 150 bps, -1 point
        If Gene Length is between 150 and 120 bps, -2 points
        If Gene Length is between 120 and 90 bps, -3 points
        If Gene Length is between 90 and 75 bps, -4 points
        If Gene Length is smaller than 75 bps, -9 points (as no gene as been recorded to be smaller than 27 amino acids long)
    '''
    score_bins = [0, 89, 119, 149, 199, 99999]
    penality = [-4, -3, -2, -1, 0]
    dataframe["length_penality"] = pd.cut(dataframe["length"], score_bins, labels=penality)


def eval_pen(dataframe):
    '''
    Create a penality system for e-values from pVOGs.
        Scoring
        -------
        If e-value is < 1E-50 bps, 3 points
        If e-value is between 1E-20 and E-50 bps, 2 point
        If e-value is between 1E-10 and 1E-20 bps, 1 points
    '''
    score_bins = [0, 1E-50, 1E-20, 1E-10, 10]
    penality = [3, 2, 1, 0]
    dataframe["e-value_penality"] = pd.cut(dataframe["e-value"], score_bins, labels=penality)


def gene_intersection(row):
    if row["strand"] == "+":
        return(pd.Index(range(row["start"], row["end"])))
    if row["strand"] == "-":
        return(pd.Index(range(row["end"], row["start"])))
    return(np.nan)


def count_overlap(dataframe):
    overlap_counter = []
    operon_counter = []
    for i in dataframe["interval"]:
        x = [len(i.intersection(y)) for y in dataframe["interval"]]
        if [1, 2, 4] in x:
            operon_counter.append(1)
        else:
            operon_counter.append(0)
        overlap_counter.append(sum(x))
    dataframe["overlap"] = abs(overlap_counter - dataframe["length"])
    dataframe["operon"] = operon_counter
    score_bins = [-999999, 10, 40, 70, 100, 999999]
    penality = [0, -1, -2, -3, -4]
    dataframe["overlap_pen"] = pd.cut(dataframe["overlap"], score_bins, labels=penality)
    dataframe["total_score"] = dataframe[["rbs_score", "score", "duplicate", "length_penality", "operon", "overlap_pen"]].fillna(0).sum(axis=1) - dataframe["truncated"].fillna(0)
    dataframe.to_csv(f"{args.output_folder}/gene_statistics.tsv", sep = "\t", index=None)
    print(f"{clr.colours.BGreen}File created at {args.output_folder}/gene_statistics.tsv{clr.colours.Colour_Off}")


def main(args):
    create_output_dir(args.output_folder)
    df = pd.read_csv(f"{args.infile}", sep="\t")
    df["length"] = df.apply(lambda row: gene_length(row), axis=1)
    length_pen(df)
    eval_pen(df)
    df["interval"] = df.apply(lambda row: gene_intersection(row), axis=1)
    count_overlap(df)


if __name__ == "__main__":
    args = my_args()
    main(args)
