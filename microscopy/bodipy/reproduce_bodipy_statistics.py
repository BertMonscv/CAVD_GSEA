#!/usr/bin/env python3
"""2026-10-12 publication companion for Figure 5b,c and Table S13A.

Input: recorded field means, not individual cell measurements.
Four field means are averaged within each of three matched experiments.
Paired two-sided t tests use n=3 (df=2); Holm adjustment covers the three
comparisons within each endpoint. Only Python's standard library is needed.
"""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics

GROUPS = ("PM", "PM+Lep", "PM+Lep+Rapa")
BATCHES = ("B1", "B2", "B3")
METRICS = ("mean_n", "mean_area")
REQUIRED = ("group", "batch", "bio_id", "field_id", "n_cells", *METRICS)


def load_fields(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(REQUIRED) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Missing columns: {', '.join(sorted(missing))}")
        rows = []
        seen = set()
        for line, row in enumerate(reader, start=2):
            if not any(row.values()):
                continue
            for key in REQUIRED:
                row[key] = row[key].strip()
            if row["group"] not in GROUPS or row["batch"] not in BATCHES:
                raise ValueError(f"Unexpected group or batch in row {line}")
            if not row["bio_id"] or not row["field_id"]:
                raise ValueError(f"Missing biological replicate or field ID in row {line}")
            key = (row["group"], row["batch"], row["field_id"])
            if key in seen:
                raise ValueError(f"Duplicate field: {key}")
            seen.add(key)
            row["n_cells"] = int(row["n_cells"])
            if row["n_cells"] <= 0:
                raise ValueError(f"Nonpositive cell count in row {line}")
            for metric in METRICS:
                row[metric] = float(row[metric])
                if not math.isfinite(row[metric]) or row[metric] < 0:
                    raise ValueError(f"Invalid {metric} in row {line}")
            rows.append(row)
    if len(rows) != 36:
        raise ValueError(f"This study has 36 fields; found {len(rows)}")
    if sum(row["n_cells"] for row in rows) != 1364:
        raise ValueError("Field cell counts do not sum to the study total of 1364")
    return rows


def paired_t_df2(first, second):
    """Exact two-sided t tail for three paired observations (df=2)."""
    difference = [a - b for a, b in zip(first, second)]
    mean = statistics.mean(difference)
    sd = statistics.stdev(difference)
    if sd == 0:
        raise ValueError("A paired test has zero variance in its differences")
    t_value = mean / (sd / math.sqrt(3))
    # Algebraically identical to 1 - |t|/sqrt(t²+2), without cancellation.
    hyp = math.hypot(t_value, math.sqrt(2))
    p_value = 2 / (hyp * (hyp + abs(t_value)))
    return t_value, p_value


def holm_adjust(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [0.0] * len(values)
    maximum = 0.0
    for rank, index in enumerate(order):
        maximum = max(maximum, values[index] * (len(values) - rank))
        adjusted[index] = min(1.0, maximum)
    return adjusted


def analyze(rows):
    experiment_values, summaries, comparisons = [], [], []
    lookup = {}
    for group, batch in itertools.product(GROUPS, BATCHES):
        subset = [r for r in rows if r["group"] == group and r["batch"] == batch]
        if len(subset) != 4 or len({r["bio_id"] for r in subset}) != 1:
            raise ValueError(f"Expected four fields from one replicate in {group}, {batch}")
        for metric in METRICS:
            value = statistics.mean(r[metric] for r in subset)
            lookup[(metric, group, batch)] = value
            experiment_values.append({"endpoint": metric, "group": group,
                                      "batch": batch, "value": value,
                                      "n_fields": 4,
                                      "n_cells": sum(r["n_cells"] for r in subset)})
    for metric in METRICS:
        for group in GROUPS:
            values = [lookup[(metric, group, batch)] for batch in BATCHES]
            summaries.append({"endpoint": metric, "group": group, "n": 3,
                              "mean": statistics.mean(values),
                              "sd": statistics.stdev(values)})
        tests = []
        for first, second in itertools.combinations(GROUPS, 2):
            t_value, p_value = paired_t_df2(
                [lookup[(metric, first, batch)] for batch in BATCHES],
                [lookup[(metric, second, batch)] for batch in BATCHES])
            tests.append({"endpoint": metric, "group_1": first,
                          "group_2": second, "n_pairs": 3, "df": 2,
                          "t": t_value, "p_raw": p_value})
        for test, adjusted in zip(tests, holm_adjust([t["p_raw"] for t in tests])):
            test["p_holm"] = adjusted
            comparisons.append(test)
    return experiment_values, summaries, comparisons


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fields", type=Path, required=True, help="Field means CSV")
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    args = parser.parse_args()
    fields = load_fields(args.fields)
    values, summaries, comparisons = analyze(fields)
    args.output.mkdir(parents=True, exist_ok=False)
    for name, rows in (("experiment_values", values), ("group_summary", summaries),
                       ("paired_holm_comparisons", comparisons)):
        write_csv(args.output / f"{name}.csv", rows)
    receipt = {"companion_date": "2026-10-12",
               "input_sha256": hashlib.sha256(args.fields.read_bytes()).hexdigest(),
               "fields": len(fields), "cells": sum(r["n_cells"] for r in fields),
               "independent_experiments_per_group": 3,
               "aggregation": "Equal mean of four field means within each experiment",
               "pairing": "Recorded B1, B2 and B3 labels across treatment groups",
               "comparison": "Two-sided paired t tests; Holm within each endpoint"}
    (args.output / "run_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"Wrote {len(values)} experiment values and {len(comparisons)} comparisons")


if __name__ == "__main__":
    main()
