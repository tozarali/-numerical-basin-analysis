"""Practical two-run diagnostic: without a reference, does the
label disagreement between two cheap runs predict the true error?"""
import numpy as np
from scipy.stats import spearmanr


def disagreement(labels_a, labels_b):
    a, b = np.asarray(labels_a), np.asarray(labels_b)
    ok = (a != 0) & (b != 0)
    return float(np.mean(a[ok] != b[ok])) if ok.any() else np.nan


def evaluate(cases, alarm=0.02, true_threshold=0.05):
    """cases: list of dicts with 'true_error' (vs reference) and
    'disagreement' (e.g. rho=1 vs rho=0.6, or h vs h/2).
    Returns rank correlation and alarm hit/false-alarm rates."""
    t = np.array([c["true_error"] for c in cases])
    d = np.array([c["disagreement"] for c in cases])
    big = t > true_threshold
    raised = d > alarm
    rho, pval = spearmanr(t, d)
    return {"spearman": float(rho), "n": int(len(t)),
            "hit_rate": float(raised[big].mean()) if big.any() else None,
            "false_alarm_rate": float(raised[~big].mean()) if (~big).any() else None,
            "alarm": alarm, "true_threshold": true_threshold,
            "note": "cases share models/ICs; p-values are descriptive only"}
