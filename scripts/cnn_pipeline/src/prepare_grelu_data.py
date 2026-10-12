#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import grelu.data.preprocess
from grelu.sequence.format import intervals_to_strings

parser = argparse.ArgumentParser(description="Extracts sequences from enhancer/background BED files and fasta assemblies for gRelu input")

parser.add_argument("--hhim_beds", required=True)
parser.add_argument("--hlat_beds", required=True)
parser.add_argument("--hhim_fasta", required=True)
parser.add_argument("--hlat_fasta", required=True)
parser.add_argument("--output", required=True)

args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

def read_bed(path):

    return pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])

def make_class(bed_path, fasta, species, label,):

    bed = read_bed(bed_path)
    
    # Extract sequences
    seqs = intervals_to_strings(bed, genome=fasta)

    df = pd.DataFrame({
        "sequence": seqs,
        "enhancer": label,
        "species": species,
    })

    return df

for split in ["train", "valid", "test"]:

    dfs = []

    for species, bed_dir, fasta in [

        ("Hhim", args.hhim_beds, args.hhim_fasta),
        ("Hlat", args.hlat_beds, args.hlat_fasta),
    ]:

        enh = make_class(f"{bed_dir}/{split}_enhancers.bed", fasta, species, 1)
        bg = make_class(f"{bed_dir}/{split}_background.bed", fasta, species, 0)

        dfs.extend([enh, bg])

    df = pd.concat(dfs, ignore_index=True)
    df = df.sample(frac=1, random_state=42).reset_index(drop=True)

    outfile = f"{args.output}/{split}.tsv"

    df.to_csv(outfile, sep="\t", index=False)

    print(f"\n{split}")
    print(df.groupby(["species", "enhancer"]).size())
    print(f"Total: {len(df)}")
