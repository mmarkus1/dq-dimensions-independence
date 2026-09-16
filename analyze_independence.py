#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = ["numpy", "scipy>=1.11"]
# ///
"""Pairwise independence analysis of data quality dimensions (Section 5.3).

Reads the problem-to-dimension assignments and computes, for every pair of
dimensions, the 2x2 contingency table over the 51 problems of the catalog,
the expected cell counts, Pearson's chi-square statistic, Fisher's exact
p-value, the Holm- and Benjamini-Hochberg-adjusted p-values, and the phi
coefficient.

Usage:
    uv run analyze_independence.py [input.csv] [output.csv]
    python3 analyze_independence.py [input.csv] [output.csv]

Requires: numpy, scipy >= 1.11 (for scipy.stats.false_discovery_control).
"""
import csv
import itertools
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_dist
from scipy.stats import false_discovery_control, fisher_exact

HERE = Path(__file__).resolve().parent
DEFAULT_IN = HERE / "problem-dimension-assignments.csv"
DEFAULT_OUT = HERE / "pairwise-independence-results.csv"


def load(path):
    """Return a list of (problem_id, set_of_affected_dimensions)."""
    problems = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            dims = set()
            primary = row["primary_affected_dimension"].strip()
            if primary:
                dims.add(primary)
            for dim in row["other_affected_dimensions"].split(","):
                dim = dim.strip()
                if dim:
                    dims.add(dim)
            problems.append((row["problem"].strip(), dims))
    return problems


def holm(pvalues):
    """Holm step-down adjusted p-values (familywise error rate)."""
    m = len(pvalues)
    adjusted = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(np.argsort(pvalues)):
        running = max(running, (m - rank) * pvalues[idx])
        adjusted[idx] = min(running, 1.0)
    return adjusted


def analyse(problems):
    n = len(problems)
    dimensions = sorted({d for _, dims in problems for d in dims})
    frequency = {d: sum(1 for _, dims in problems if d in dims) for d in dimensions}

    rows = []
    for a, b in itertools.combinations(dimensions, 2):
        n11 = sum(1 for _, dims in problems if a in dims and b in dims)
        n10 = sum(1 for _, dims in problems if a in dims and b not in dims)
        n01 = sum(1 for _, dims in problems if a not in dims and b in dims)
        n00 = n - n11 - n10 - n01

        row_a, row_b = n11 + n10, n01 + n00
        col_a, col_b = n11 + n01, n10 + n00
        expected = [
            [row_a * col_a / n, row_a * col_b / n],
            [row_b * col_a / n, row_b * col_b / n],
        ]
        observed = [[n11, n10], [n01, n00]]

        chi2 = sum(
            (observed[i][j] - expected[i][j]) ** 2 / expected[i][j]
            for i in range(2)
            for j in range(2)
            if expected[i][j] > 0
        )
        _, fisher_p = fisher_exact(observed)
        denominator = row_a * row_b * col_a * col_b
        phi = (n11 * n00 - n10 * n01) / math.sqrt(denominator) if denominator else float("nan")

        rows.append(
            {
                "dimension_a": a,
                "dimension_b": b,
                "freq_a": frequency[a],
                "freq_b": frequency[b],
                "n11": n11,
                "n10": n10,
                "n01": n01,
                "n00": n00,
                "exp_11": expected[0][0],
                "exp_10": expected[0][1],
                "exp_01": expected[1][0],
                "exp_00": expected[1][1],
                "min_expected": min(min(r) for r in expected),
                "chi2": chi2,
                "chi2_p": chi2_dist.sf(chi2, 1),
                "fisher_p": fisher_p,
                "phi": phi,
            }
        )

    pvalues = np.array([r["fisher_p"] for r in rows])
    for row, h, bh in zip(rows, holm(pvalues), false_discovery_control(pvalues, method="bh")):
        row["fisher_p_holm"] = h
        row["fisher_p_bh"] = bh

    return rows, frequency, n


def write(rows, path):
    columns = [
        "dimension_a", "dimension_b", "freq_a", "freq_b",
        "n11", "n10", "n01", "n00",
        "exp_11", "exp_10", "exp_01", "exp_00", "min_expected",
        "chi2", "chi2_p", "fisher_p", "fisher_p_holm", "fisher_p_bh", "phi",
    ]
    integers = {"freq_a", "freq_b", "n11", "n10", "n01", "n00"}
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        for row in sorted(rows, key=lambda r: -abs(r["phi"])):
            writer.writerow(
                [row[c] if c in integers or isinstance(row[c], str) else f"{row[c]:.6f}"
                 for c in columns]
            )


def sensitivity(rows, threshold=5.0):
    """Compare using Fisher throughout against a mixed strategy.

    The mixed strategy applies Pearson's chi-square test wherever its
    expected-count assumption holds and Fisher's exact test elsewhere. Both
    families of p-values are Holm-corrected over all comparisons.
    """
    mixed = np.array([
        r["chi2_p"] if r["min_expected"] >= threshold else r["fisher_p"] for r in rows
    ])
    uniform = np.array([r["fisher_p"] for r in rows])
    mixed_holm, uniform_holm = holm(mixed), holm(uniform)

    print("\nsensitivity to the choice of test")
    print(f"  Fisher throughout, uncorrected p < 0.05: {int((uniform < 0.05).sum())}")
    print(f"  mixed strategy,    uncorrected p < 0.05: {int((mixed < 0.05).sum())}")
    print(f"  Fisher throughout, Holm-adjusted p < 0.05: {int((uniform_holm < 0.05).sum())}")
    print(f"  mixed strategy,    Holm-adjusted p < 0.05: {int((mixed_holm < 0.05).sum())}")
    same = bool(((uniform_holm < 0.05) == (mixed_holm < 0.05)).all())
    print(f"  identical set of significant pairs after correction: {same}")

    admissible = [r for r in rows if r["min_expected"] >= threshold]
    deltas = [r["fisher_p"] - r["chi2_p"] for r in admissible]
    print(f"  pairs meeting the expected-count criterion: {len(admissible)}")
    print(f"  Fisher is more conservative on all of them: {all(d > 0 for d in deltas)}")
    print(f"  mean difference in p: {sum(deltas)/len(deltas):+.4f}")


def main():
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_IN
    target = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_OUT

    problems = load(source)
    rows, frequency, n = analyse(problems)

    print(f"problems: {n}")
    print(f"dimensions with at least one assignment: {len(frequency)}")
    print(f"assignments: {sum(frequency.values())}")
    print(f"pairs: {len(rows)}")
    print(f"pairs with min expected count < 5: {sum(1 for r in rows if r['min_expected'] < 5)}")
    print(f"pairs with min expected count < 1: {sum(1 for r in rows if r['min_expected'] < 1)}")
    print(f"uncorrected Fisher p < 0.05: {sum(1 for r in rows if r['fisher_p'] < 0.05)}")
    print(f"Holm-adjusted p < 0.05: {sum(1 for r in rows if r['fisher_p_holm'] < 0.05)}")
    print(f"Benjamini-Hochberg q < 0.05: {sum(1 for r in rows if r['fisher_p_bh'] < 0.05)}")

    sensitivity(rows)

    phis = np.array([abs(r["phi"]) for r in rows])
    print(f"median |phi|: {np.median(phis):.3f}")
    print(f"pairs with |phi| < 0.3: {int((phis < 0.3).sum())}")

    write(rows, target)
    print(f"written: {target}")


if __name__ == "__main__":
    main()
