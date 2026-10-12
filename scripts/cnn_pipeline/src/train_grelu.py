#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import shutil
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
)
from grelu.data.dataset import DFSeqDataset
from grelu.lightning import LightningModel

parser = argparse.ArgumentParser(description="Train gReLU enhancer classifier")
parser.add_argument("--train", required=True, help="Training TSV")
parser.add_argument("--valid", required=True, help="Validation TSV")
parser.add_argument("--test", required=True, help="Test TSV")
parser.add_argument("--output", required=True, help="Output directory")

args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

# load sequence tables

print("\nLoading data")

train = pd.read_csv(args.train, sep="\t")
valid = pd.read_csv(args.valid, sep="\t")
test = pd.read_csv(args.test, sep="\t")

print(f"Train: {len(train)}")
print(f"Valid: {len(valid)}")
print(f"Test:  {len(test)}")

print("\nTraining class counts:")
print(train.groupby(["species", "enhancer"]).size())

print("\nValidation class counts:")
print(valid.groupby(["species", "enhancer"]).size())

print("\nTest class counts:")
print(test.groupby(["species", "enhancer"]).size())

# cross-species comparison

train = (train[train["species"] == "Hlat"].reset_index(drop=True))
valid = (valid[valid["species"] == "Hhim"].reset_index(drop=True))
 test = (test[ test["species"] == "Hhim"].reset_index(drop=True)) 
 
print("\nCross-species dataset:")
 
print(
        f"Train Hlat: {len(train)} "
        f"({train['enhancer'].sum()} enhancers)"
)
 
print(
        f"Valid Hhim: {len(valid)} "
        f"({valid['enhancer'].sum()} enhancers)"
)
 
print(
        f"Test Hhim: {len(test)} "
        f"({test['enhancer'].sum()} enhancers)"
)
 
assert set(train["species"]) == {"Hlat"}
assert set(valid["species"]) == {"Hhim"}
assert set(test["species"]) == {"Hhim"}
assert len(train) > 0
assert len(valid) > 0
assert len(test) > 0
 
# gReLU datasets

train_ds = DFSeqDataset(train[["sequence", "enhancer"]], seq_len=501, rc=True, max_seq_shift=0)
valid_ds = DFSeqDataset(valid[["sequence", "enhancer"]], seq_len=501)
test_ds = DFSeqDataset(test[["sequence", "enhancer"]], seq_len=501)

print("\ngReLU datasets:")
print(f"Train: {len(train_ds)}")
print(f"Valid: {len(valid_ds)}")
print(f"Test:  {len(test_ds)}")

# model

model_params = {
    "model_type": "DilatedConvModel", "n_tasks": 1, "channels": 128, "stem_kernel_size": 21, "kernel_size": 5,
    "dilation_mult": 2, "n_conv": 8, "crop_len": 0, "final_pool_func": "avg"
}

# training parameters

train_params = {
    "task": "binary", "lr": 8.609554010119795e-05, "batch_size": 256, "max_epochs": 100, "optimizer": "adam",
    "logger": "csv", "checkpoint": True, "devices": 0, "num_workers": 8, "save_dir": args.output
}

# construct model

print("\nBuilding gReLU model")

model = LightningModel(model_params=model_params, train_params=train_params,)

# train

print("\nTraining")

trainer = model.train_on_dataset(train_dataset=train_ds, val_dataset=valid_ds,)

# load best checkpoint

best_checkpoint = trainer.checkpoint_callback.best_model_path

if not best_checkpoint:
    raise RuntimeError("No checkpoint was saved")

print(f"\nBest checkpoint:\n{best_checkpoint}")

best_model_path = os.path.join(args.output, "best.ckpt")
shutil.copy2(best_checkpoint, best_model_path)

print(f"Saved best model: {best_model_path}")


# load trained model
model = LightningModel.load_from_checkpoint(best_checkpoint, weights_only=False)

# prediction helper

def get_predictions(model, dataset):

    pred = model.predict_on_dataset(dataset, devices=0, num_workers=8, batch_size=256)
    pred = np.asarray(pred).reshape(-1)

    return pred

# validation predictions

print("\nPredicting validation set")

valid_scores = get_predictions(model, valid_ds)

valid_y = valid["enhancer"].to_numpy()

# determine threshold from validation

TARGET_PRECISION = 0.90

pr_precision, pr_recall, pr_thresholds = (precision_recall_curve(valid_y, valid_scores))

# precision/recal arrays have one extra element compared with thresholds

eligible = np.where(pr_precision[:-1] >= TARGET_PRECISION)[0]

if len(eligible) == 0:
    raise RuntimeError(
        f"No validation threshold achieves "
        f"{TARGET_PRECISION:.2f} precision"
    )

# thresholds satisfying target precision choise highest recall

best_idx = eligible[np.argmax(pr_recall[:-1][eligible])]
threshold = float(pr_thresholds[best_idx])

print("\nValidation threshold selection")
print(f"Target precision: {TARGET_PRECISION:.4f}")
print(f"Threshold:        {threshold:.6f}")
print(
    f"Precision:        "
    f"{pr_precision[best_idx]:.4f}"
)
print(
    f"Recall:           "
    f"{pr_recall[best_idx]:.4f}"
)

# test predictions

print("\nPredicting held-out test set")

test_scores = get_predictions(model, test_ds)

test_y = test["enhancer"].to_numpy()

# evaluation function

def evaluate_split(split, y_true, scores, threshold):

    y_pred = (scores >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    specificity = (tn / (tn + fp)
        if (tn + fp) > 0
        else np.nan
    )

    npv = (tn / (tn + fn)
        if (tn + fn) > 0
        else np.nan
    )

    return {
        "split": split, "n_sequences": len(y_true), "n_enhancers": int(np.sum(y_true == 1)), "n_background": int(np.sum(y_true == 0)),
        "auroc": roc_auc_score(y_true, scores), "auprc": average_precision_score(y_true, scores), "threshold": threshold,
        "accuracy": accuracy_score(y_true, y_pred), "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0), "f1": f1_score(y_true, y_pred, zero_division=0),
        "specificity": specificity, "npv": npv, "tp": int(tp), "fp": int(fp), "tn": int(tn), "fn": int(fn),
        "mean_enhancer_score": float(scores[y_true == 1].mean()), "mean_background_score": float(scores[y_true == 0].mean())}

# evaluate valid and test 

valid_metrics = evaluate_split("valid", valid_y, valid_scores, threshold)
test_metrics = evaluate_split("test", test_y, test_scores, threshold)

# save a summary TSV

summary = pd.DataFrame([valid_metrics, test_metrics,])
summary_path = os.path.join(args.output, "evaluation_summary.tsv")
summary.to_csv(summary_path, sep="\t", index=False, float_format="%.6f")

# save per sequence validation predictions

valid_output = valid.copy()

valid_output["enhancer_score"]= valid_scores
valid_output["predicted_class"] = (valid_scores >= threshold).astype(int)
valid_output.to_csv(os.path.join(args.output, "validation_predictions.tsv"), sep="\t", index=False, float_format="%.6f")

# save per sequence test predictions

test_output = test.copy()

test_output["enhancer_score"] = test_scores
test_output["predicted_class"] = (test_scores >= threshold).astype(int)
test_output.to_csv(os.path.join(args.output, "test_predictions.tsv"), sep="\t", index=False, float_format="%.6f")

# print final summary

print("\n============================================================")
print("Final model performance")
print("============================================================")

print(summary.to_string(index=False))

print(f"\nThreshold determined from validation: {threshold:.6f}")
print(f"Summary: {summary_path}")
print(f"Model:   {best_model_path}")


