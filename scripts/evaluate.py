import pandas as pd
import numpy as np
from sklearn.metrics import matthews_corrcoef, confusion_matrix

def parse_tbl(path):
    """Parse hmmsearch --tblout. Returns dict {uniprot_acc: best_evalue}."""
    hits = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 5: continue
            target = parts[0]                  # sp|P12345|NAME_SPECIES
            acc = target.split("|")[1] if "|" in target else target
            evalue = float(parts[4])
            if acc not in hits or evalue < hits[acc]:
                hits[acc] = evalue
    return hits

def read_ids(fasta):
    ids = []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split("|")[1])
    return ids

def evaluate(fold):
    pos_ids = read_ids(f"/workspace/validation/positives_set{fold}.fasta")
    neg_ids = read_ids(f"/workspace/validation/negatives_set{fold}.fasta")
    pos_hits = parse_tbl(f"positives_set{fold}.tbl")
    neg_hits = parse_tbl(f"negatives_set{fold}.tbl")

    # build score table: every id, its best E-value (1e6 = no hit)
    rows = []
    for aid in pos_ids:
        rows.append((aid, 1, pos_hits.get(aid, 1e6)))
    for aid in neg_ids:
        rows.append((aid, 0, neg_hits.get(aid, 1e6)))
    df = pd.DataFrame(rows, columns=["id","label","evalue"])
    return df

def scan_thresholds(df):
    results = []
    # full-sequence E-value thresholds to test
    thresholds = [10**i for i in range(-30, 5)]
    for t in thresholds:
        pred = (df.evalue <= t).astype(int)
        tn, fp, fn, tp = confusion_matrix(df.label, pred, labels=[0,1]).ravel()
        mcc = matthews_corrcoef(df.label, pred) if (tp+fp>0 and tn+fn>0) else 0
        acc = (tp+tn)/(tp+tn+fp+fn)
        prec = tp/(tp+fp) if tp+fp else 0
        rec = tp/(tp+fn) if tp+fn else 0
        results.append((t, tp, fp, fn, tn, acc, prec, rec, mcc))
    return pd.DataFrame(results, columns=["E_threshold","TP","FP","FN","TN","Acc","Prec","Recall","MCC"])

# 2-fold CV
df1 = evaluate(1)
df2 = evaluate(2)

print("=== FOLD 1: train on set1 metrics, test on set2 ===")
# Use set1 to find best threshold, evaluate on set2
scan1 = scan_thresholds(df1)
best_t1 = scan1.loc[scan1.MCC.idxmax(), "E_threshold"]
print(f"Best threshold from set1: {best_t1:.1e}, MCC={scan1.MCC.max():.4f}")

print("\n=== FOLD 2: train on set2, test on set1 ===")
scan2 = scan_thresholds(df2)
best_t2 = scan2.loc[scan2.MCC.idxmax(), "E_threshold"]
print(f"Best threshold from set2: {best_t2:.1e}, MCC={scan2.MCC.max():.4f}")

# Apply set1-derived threshold to set2 (and vice versa) — true CV
def report(df, t, name):
    pred = (df.evalue <= t).astype(int)
    tn, fp, fn, tp = confusion_matrix(df.label, pred, labels=[0,1]).ravel()
    mcc = matthews_corrcoef(df.label, pred)
    acc = (tp+tn)/(tp+tn+fp+fn)
    prec = tp/(tp+fp) if tp+fp else 0
    rec = tp/(tp+fn) if tp+fn else 0
    print(f"{name}: TP={tp} FP={fp} FN={fn} TN={tn} Acc={acc:.5f} Prec={prec:.4f} Rec={rec:.4f} MCC={mcc:.4f}")

print("\n=== 2-fold CV results ===")
report(df2, best_t1, "Set2 evaluated with set1 threshold")
report(df1, best_t2, "Set1 evaluated with set2 threshold")

# combined
df_all = pd.concat([df1, df2])
print("\n=== Threshold scan on all data ===")
scan_all = scan_thresholds(df_all)
print(scan_all.to_string(index=False))
scan_all.to_csv("threshold_scan.csv", index=False)

# also save full score table
df_all.to_csv("all_scores.csv", index=False)
print("\nSaved: threshold_scan.csv, all_scores.csv")
