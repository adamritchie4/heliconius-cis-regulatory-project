#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd
import pysam
import torch
from grelu.lightning import LightningModel

parser = argparse.ArgumentParser()
parser.add_argument("--fasta", required=True)
parser.add_argument("--bed", required=True)
parser.add_argument("--model", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--threshold", type=float, default=0.804347)
parser.add_argument("--batch_size", type=int, default=4096)
args = parser.parse_args()

# load model

print(f"Loading model: {args.model}")

model = LightningModel.load_from_checkpoint(args.model, weights_only=False)
model.eval()

# open genome

genome = pysam.FastaFile(args.fasta)

seqs = []
coords = []

total = 0
positive = 0

# prediction function

def predict_batch(seqs, coords, out):

    global positive
    scores = np.asarray(model.predict_on_seqs(seqs, device=0)).reshape(-1)

    for (chrom, start, end), score in zip(coords, scores):

        if score >= args.threshold:

            out.write(f"{chrom}\t" f"{start}\t" f"{end}\t" f"{score:.6f}\n")
            positive += 1

# read BED and predict

with open(args.bed) as bed, open(args.output, "w") as out:

    for line in bed:

        chrom, start, end = line.rstrip().split()[:3]
        start = int(start)
        end = int(end)
        seq = genome.fetch(chrom, start, end).upper()

        seqs.append(seq)
        coords.append((chrom, start, end))

        total += 1

        # submit in batch

        if len(seqs) >= args.batch_size:

            predict_batch(seqs, coords, out)

            seqs = []
            coords = []

            if total % 100000 == 0:
                print(
                    f"Processed: {total:,} "
                    f"Positive: {positive:,}"
                )

    # finish any remaining partial batch

    if seqs:

        predict_batch(seqs, coords, out)

genome.close()

print("\nPrediction complete")
print(f"Sequences processed: {total:,}")
print(f"Positive predictions: {positive:,}")
print(f"Threshold: {args.threshold}")
print(f"Output: {args.output}")
