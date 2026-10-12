#!/usr/bin/env python3

import argparse
import os
import re
import yaml
import pandas as pd

parser = argparse.ArgumentParser(description=("Split enhancer/background BED files by chromosome while retaining original assembly sequence names."))

parser.add_argument("--species", required=True, help="Species abbreviation, e.g. Hhim or Hlat")
parser.add_argument("--config", required=True, help="YAML file containing train/valid/test chromosome splits")
parser.add_argument("--input", required=True, help="Directory containing enhancers.bed and background.bed")
parser.add_argument("--output", required=True, help="Output directory")

args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

# convert assembly partition names to parent chromosomes

def get_split_chrom(chrom):

    if args.species == "Hlat":

        # Example:
        # Hel_chr11_1 -> Hlat11
        # Hel_chr11_2 -> Hlat11

        match = re.fullmatch(
            r"Hel_chr(\d+)_\d+",
            chrom
        )

        if match:
            chromosome = int(match.group(1))
            return f"Hlat{chromosome:02d}"

        if re.fullmatch(r"Hlat\d+", chrom):
            return chrom

    # Hhim chromosome names already match split names
    return chrom

# load chromosome splits

with open(args.config) as f:
    splits = yaml.safe_load(f)


if args.species not in splits:
    raise ValueError(
        f"Species '{args.species}' not found in {args.config}"
    )

train = set(splits[args.species]["train"])
valid = set(splits[args.species]["valid"])
test = set(splits[args.species]["test"])

# ensure chromosomes do not occur in multiple sets

if train & valid:
    raise ValueError(
        f"Chromosomes occur in both train and valid: "
        f"{sorted(train & valid)}"
    )

if train & test:
    raise ValueError(
        f"Chromosomes occur in both train and test: "
        f"{sorted(train & test)}"
    )

if valid & test:
    raise ValueError(
        f"Chromosomes occur in both valid and test: "
        f"{sorted(valid & test)}"
    )

print(f"Species: {args.species}")
print(f"Train chromosomes: {sorted(train)}")
print(f"Valid chromosomes: {sorted(valid)}")
print(f"Test chromosomes:  {sorted(test)}")

# process enhancer and background classes

for dataset in ["enhancers", "background"]:

    print(f"\nProcessing {dataset}")

    infile = os.path.join(args.input, f"{dataset}.bed")

    if not os.path.isfile(infile):
        raise FileNotFoundError(f"Could not find input BED: {infile}")

    df = pd.read_csv(infile, sep="\t", header=None, names=["chrom", "start", "end"], usecols=[0, 1, 2])

    # temp chromosomes for splitting

    df["split_chrom"] = df["chrom"].apply(get_split_chrom)

    # split by parent chromosome

    train_df = df[df["split_chrom"].isin(train)].copy()
    valid_df = df[df["split_chrom"].isin(valid)].copy()
    test_df = df[ df["split_chrom"].isin(test) ].copy()

    assigned = (df["split_chrom"].isin(train) | df["split_chrom"].isin(valid) | df["split_chrom"].isin(test))

    unassigned = df[~assigned]

    if len(unassigned) > 0:
        print(
                f"WARNING: {len(unassigned)} {dataset} regions "
                "were not assigned to any split."
        )
        print(
            "Unassigned parent chromosomes:",
            sorted(
                unassigned["split_chrom"].unique()
            )
        )

        for split_name, split_df in [
        ("train", train_df),
        ("valid", valid_df),
        ("test", test_df)
    ]:

        outfile = os.path.join(args.output, f"{split_name}_{dataset}.bed")
        split_df[["chrom", "start", "end"]].to_csv(outfile, sep="\t", header=False, index=False)

    print(
        f"Total: {len(df)} | "
        f"Train: {len(train_df)} | "
        f"Valid: {len(valid_df)} | "
        f"Test: {len(test_df)} | "
        f"Unassigned: {len(unassigned)}"
    )

print("\nFinished splitting BED files.")
