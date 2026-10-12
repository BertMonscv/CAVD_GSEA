#!/usr/bin/env python3
"""2026-10-12 publication companion: Figure 5e-h/Table S13B-C statistics.

Read recorded experiment-level values. This does not recalculate cell ROIs,
thresholds, Manders coefficients, or per-cell enrichment from source images.
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
EXPERIMENTS = ("1", "2", "3")
PANELS = {"e": "fracG", "f": "fracR", "g": "M2", "h": "enrich"}
REQUIRED = ("panel", "endpoint", "group", "matched_experiment", "value", "unit")


def load_values(path):
    result = {}
    units = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(REQUIRED) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Missing columns: {sorted(missing)}")
        for row in reader:
            if row["panel"] not in PANELS:
                continue  # The complete Figure 5 table also contains panels b/c.
            if row["endpoint"] != PANELS[row["panel"]]:
                raise ValueError("Panel and endpoint labels do not match")
            group, experiment = row["group"], row["matched_experiment"]
            if group not in GROUPS or experiment not in EXPERIMENTS:
                raise ValueError("Unknown group or matched experiment label")
            key = (row["endpoint"], group, experiment)
            if key in result:
                raise ValueError(f"Duplicate experiment value: {key}")
            value = float(row["value"])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid measurement: {key}")
            result[key] = value
            endpoint = row["endpoint"]
            units.setdefault(endpoint, row["unit"])
            if units[endpoint] != row["unit"]:
                raise ValueError(f"Mixed units within {endpoint}")
    expected = set(itertools.product(PANELS.values(), GROUPS, EXPERIMENTS))
    if set(result) != expected:
        raise ValueError(f"Expected 36 experiment values; missing {sorted(expected-set(result))}")
    expected_units = {"fracG": "percent", "fracR": "percent",
                      "M2": "dimensionless fraction", "enrich": "dimensionless ratio"}
    if units != expected_units:
        raise ValueError(f"Unexpected units: {units}")
    return result, units


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(values, units):
    summaries, comparisons = [], []
    for panel, endpoint in PANELS.items():
        for group in GROUPS:
            sample = [values[(endpoint, group, e)] for e in EXPERIMENTS]
            summaries.append({"panel": panel, "endpoint": endpoint, "group": group,
                              "n": 3, "mean": statistics.mean(sample),
                              "sd": statistics.stdev(sample), "unit": units[endpoint]})
        for first, second in itertools.combinations(GROUPS, 2):
            differences = [values[(endpoint, first, e)] - values[(endpoint, second, e)]
                           for e in EXPERIMENTS]
            sd = statistics.stdev(differences)
            if sd == 0:
                raise ValueError("A paired test has zero variance in its differences")
            t_value = statistics.mean(differences) / (sd / math.sqrt(3))
            hyp = math.hypot(t_value, math.sqrt(2))
            p_value = 2 / (hyp * (hyp + abs(t_value)))
            comparisons.append({"panel": panel, "endpoint": endpoint,
                                "group_1": first, "group_2": second,
                                "n_pairs": 3, "df": 2, "t": t_value,
                                "p_unadjusted": p_value,
                                "comparison_status": "exploratory" if
                                (first, second) == ("PM", "PM+Lep+Rapa") else "pre-specified"})
    return summaries, comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--values", required=True, type=Path,
                        help="Recorded experiment-level Figure 5 values CSV")
    parser.add_argument("--output", required=True, type=Path, help="New output directory")
    args = parser.parse_args()
    values, units = load_values(args.values)
    summaries, comparisons = analyze(values, units)
    args.output.mkdir(parents=True, exist_ok=False)
    write_csv(args.output / "group_summary.csv", summaries)
    write_csv(args.output / "paired_comparisons.csv", comparisons)
    receipt = {"companion_date": "2026-10-12", "input_level": "experiment means",
               "input_sha256": hashlib.sha256(args.values.read_bytes()).hexdigest(),
               "experiment_values": len(values), "n_per_group": 3,
               "comparison": "Two-sided paired t tests, df=2, unadjusted P values",
               "enrichment": "Uses recorded means of per-cell ratios, not a ratio of experiment means"}
    (args.output / "run_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"Wrote {len(summaries)} group summaries and {len(comparisons)} comparisons")


if __name__ == "__main__":
    main()
