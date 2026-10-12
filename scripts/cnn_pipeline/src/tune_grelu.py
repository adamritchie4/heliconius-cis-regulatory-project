#!/usr/bin/env python3

import os
import gc
import argparse
import pandas as pd
import torch
import optuna
import grelu.lightning
from grelu.data.dataset import DFSeqDataset

parser = argparse.ArgumentParser("Uses Optuna for determining best hyperparameters for DilatedConvModel training")

parser.add_argument("--train", required=True)
parser.add_argument("--valid", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--trials", type=int, default=20)

args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

# load data

print("\nLoading data")

train = pd.read_csv(args.train, sep="\t")
valid = pd.read_csv(args.valid, sep="\t")

print(f"Train: {len(train)}")
print(f"Valid: {len(valid)}")

# gReLU datasets

train_ds = DFSeqDataset(train[["sequence", "enhancer"]], seq_len=500, rc=True, max_seq_shift=0, augment_mode="serial")
valid_ds = DFSeqDataset(valid[["sequence", "enhancer"]], seq_len=500)

# Optuna objective

def objective(trial):

    print("\n" + "=" * 50)
    print(f"Trial {trial.number}")
    print("=" * 50)

    # parameters to tune

    stem_kernel_size = trial.suggest_categorical("stem_kernel_size", [15,21,31])
    kernel_size = trial.suggest_categorical("kernel_size", [3,5])
    dilation_mult = trial.suggest_categorical("dilation_mult", [1,2])
    lr = trial.suggest_float("lr", 3e-5, 3e-4, log=True)
    batch_size = trial.suggest_categorical("batch_size", [128,256,512])
    final_pool_func = trial.suggest_categorical("final_pool_func", ["avg","max"])

    # fixed architecture + trial parameters

    model_params = {

        "model_type": "DilatedConvModel",
        "n_tasks": 1,

        # parameters
        "channels": 128,
        "stem_kernel_size": stem_kernel_size,
        "kernel_size": kernel_size,
        "dilation_mult": dilation_mult,
        "n_conv": 8,
        "crop_len": 0,
        "final_pool_func": final_pool_func
    }

    trial_dir = os.path.join(args.output, f"trial_{trial.number:03d}")
    os.makedirs(trial_dir, exist_ok=True)
    train_params = {
            "task": "binary", "lr": lr, "batch_size": batch_size, "max_epochs": 70, "optimizer": "adam",
            "logger": "csv", "checkpoint": True, "devices": 0, "num_workers": 8, "save_dir": trial_dir
    }


    print("\nModel parameters:")
    print(model_params)

    print("\nTraining parameters:")
    print(train_params)


    # build and train

    model = grelu.lightning.LightningModel(model_params=model_params, train_params=train_params)
    trainer = model.train_on_dataset(train_dataset=train_ds, val_dataset=valid_ds,)

    # read validation metrics

    metrics_file = os.path.join(trainer.logger.log_dir, "metrics.csv")
    metrics = pd.read_csv(metrics_file)

    val_auprc = (metrics["val_avgprec"].dropna())

    if len(val_auprc) == 0:
        raise RuntimeError("No val_avgprec values found")

    best_val_auprc = float(val_auprc.max())
    val_auroc = (metrics["val_auroc"].dropna())
    best_val_auroc = float(val_auroc.max())

    print(
        f"\nBest validation AUPRC: "
        f"{best_val_auprc:.6f}"
    )

    print(
        f"Best validation AUROC: "
        f"{best_val_auroc:.6f}"
    )

    trial.set_user_attr("best_val_auroc", best_val_auroc)
    return best_val_auprc

# run search
sampler = optuna.samplers.TPESampler(seed=42)

storage = (f"sqlite:///" f"{os.path.abspath(args.output)}" f"/optuna.db")     
study = optuna.create_study(study_name="heliconius_grelu_500bp", direction="maximize", sampler=sampler, storage=storage, load_if_exists=True)
study.enqueue_trial({"stem_kernel_size": 21, "kernel_size": 3, "dilation_mult": 2, "lr": 1e-4, "batch_size": 256, "final_pool_func": "avg"})
study.optimize(objective, n_trials=args.trials, gc_after_trial=True,)

# Save results

print("\n" + "=" * 60)
print("Hyperparameter search complete")
print("=" * 60)

print(
    f"\nBest validation AUPRC: "
    f"{study.best_value:.6f}"
)

print("\nBest parameters:")

for key, value in (
    study.best_params.items()
):

    print(
        f"{key}: {value}"
    )


results = study.trials_dataframe()
results.to_csv(os.path.join(args.output, "hyperparameter_trials.tsv"), sep="\t", index=False)
