import argparse
import re
import pandas as pd

parser = argparse.ArgumentParser(description="Generate a 1:1 genomic background matched only for N content.")

parser.add_argument("--enhancers", required=True, help="bedtools nuc output for enhancer sequences")
parser.add_argument("--background", required=True, help="bedtools nuc output for background candidate sequences")
parser.add_argument("--output", required=True, help="Output BED file")
parser.add_argument("--max_n_fraction", type=float, default=0.01, help="Maximum permitted N fraction")
parser.add_argument("--seed", type=int, default=42, help="Random seed")

args = parser.parse_args()

def load_nuc_file(filepath):
    with open(filepath) as f:
        header = f.readline().strip().split("\t")

    header = [
        re.sub(r"^[0-9]+_", "", col.lstrip("#"))
        for col in header
    ]

    df = pd.read_csv(filepath, sep="\t", comment="#", header=None)

    columns = {
        name: i
        for i, name in enumerate(header)
    }

    for col in ["num_N", "seq_len"]:
        if col not in columns:
            raise ValueError(
                f"Missing column '{col}' in {filepath}\n"
                f"Available columns: {header}"
            )

    return df, columns

enh, enh_cols = load_nuc_file(args.enhancers)
bg, bg_cols = load_nuc_file(args.background)

enh_n = enh_cols["num_N"]
enh_len = enh_cols["seq_len"]
bg_n = bg_cols["num_N"]
bg_len = bg_cols["seq_len"]

enh[enh_n] = pd.to_numeric(enh[enh_n])
enh[enh_len] = pd.to_numeric(enh[enh_len])
bg[bg_n] = pd.to_numeric(bg[bg_n])
bg[bg_len] = pd.to_numeric(bg[bg_len])

enh["N_fraction"] = enh[enh_n] / enh[enh_len]
bg["N_fraction"] = bg[bg_n] / bg[bg_len]

enh = enh[enh["N_fraction"] <= args.max_n_fraction].copy()
bg = bg[bg["N_fraction"] <= args.max_n_fraction].copy()

print(f"Enhancers after N filtering: {len(enh):,}")
print(f"Background candidates after N filtering: {len(bg):,}")

matched = []

for n_count, enh_group in enh.groupby(enh_n):
    candidates = bg[bg[bg_n] == n_count]

    if len(candidates) == 0:
        print(f"Warning: no background candidates with num_N={n_count}")
        continue

    n_required = min(len(enh_group), len(candidates))

    if n_required < len(enh_group):
        print(
            f"Warning: num_N={n_count}: needed {len(enh_group):,}, "
            f"available {len(candidates):,}"
        )

    matched.append(candidates.sample(n=n_required, random_state=args.seed))

if not matched:
    raise RuntimeError("No N-matched background sequences found.")

background = pd.concat(matched, ignore_index=True)

print()
print(f"Enhancers : {len(enh):,}")
print(f"Background: {len(background):,}")
print(f"Mean enhancer N fraction : {enh['N_fraction'].mean():.6f}")
print(f"Mean background N fraction: {background['N_fraction'].mean():.6f}")

print("\nEnhancer N counts:")
print(enh[enh_n].value_counts().sort_index().to_string())

print("\nBackground N counts:")
print(background[bg_n].value_counts().sort_index().to_string())

background.iloc[:, 0:3].to_csv(args.output, sep="\t", header=False, index=False)

print(f"\nWritten {args.output}")
