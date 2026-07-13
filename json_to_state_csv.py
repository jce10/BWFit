#!/usr/bin/env python3
"""
Extract state fit results from BWFit JSON summary files and write a tidy CSV.

Output columns match the earlier summary CSV:
permutation,fit,theta,state,shape,yield_set,yield_fit,yield_err,
centroid_set,centroid_fit,centroid_err,width_set,width_fit,width_err,
phase,phase_err,chi2,ndf,chi2_ndf

Examples
--------
# One directory of JSON files
python json_to_state_csv.py results/json_fits -o free_centroid_free_width_all_states.csv

# Explicit files and a custom permutation label
python json_to_state_csv.py 10deg.json 15deg.json 20deg.json \
    -o free_centroid_free_width_all_states.csv \
    --permutation free_centroid_free_width
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any, Iterable

COLUMNS = [
    "permutation",
    "fit",
    "theta",
    "state",
    "shape",
    "yield_set",
    "yield_fit",
    "yield_err",
    "centroid_set",
    "centroid_fit",
    "centroid_err",
    "width_set",
    "width_fit",
    "width_err",
    "phase",
    "phase_err",
    "chi2",
    "ndf",
    "chi2_ndf",
]


def as_number(value: Any, default: float | int | None = None) -> Any:
    """Return a JSON value as a number when possible."""
    if value is None:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def get_nested(data: dict[str, Any], *keys: str, default: Any = None) -> Any:
    """Safely walk nested dictionaries."""
    current: Any = data
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def normalize_param(param: dict[str, Any] | None) -> tuple[Any, Any, Any]:
    """Return set, fit, err from a BWFit JSON parameter block."""
    if not isinstance(param, dict):
        return None, None, None
    return (
        as_number(param.get("set")),
        as_number(param.get("fit")),
        as_number(param.get("err")),
    )


def theta_from_filename(path: Path) -> float | None:
    """Try to recover theta from filenames like 10deg_fit.json or theta_15.json."""
    match = re.search(r"(?:^|[_-])(\d+(?:\.\d+)?)(?:deg|degree|degrees)(?:[_-]|$)", path.stem, re.I)
    if match:
        return float(match.group(1))

    match = re.search(r"theta[_-]?(\d+(?:\.\d+)?)", path.stem, re.I)
    if match:
        return float(match.group(1))

    return None


def json_paths(inputs: Iterable[Path], recursive: bool) -> list[Path]:
    """Expand input files/directories into a sorted unique list of JSON files."""
    paths: list[Path] = []
    for item in inputs:
        if item.is_dir():
            globber = item.rglob if recursive else item.glob
            paths.extend(globber("*.json"))
        elif item.is_file() and item.suffix.lower() == ".json":
            paths.append(item)
        else:
            print(f"[skip] Not a JSON file or directory: {item}")

    return sorted(set(paths), key=lambda p: str(p))


def extract_rows(path: Path, permutation: str, fit_number: int) -> list[dict[str, Any]]:
    """Extract one row per state from one BWFit JSON summary."""
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    fit_config = data.get("fit_config", {})
    fit_result = data.get("fit_result", {})

    theta = as_number(fit_config.get("angle_deg"), theta_from_filename(path))

    quality = fit_result.get("quality", {})
    chi2 = as_number(quality.get("chi2"))
    ndf_raw = quality.get("ndf")
    ndf = int(ndf_raw) if ndf_raw is not None else None
    chi2_ndf = as_number(quality.get("chi2_ndf"))

    interference = fit_result.get("interference", {})
    phase = as_number(interference.get("phase_fit_rad"))
    phase_err = as_number(interference.get("phase_err_rad"))

    states = fit_result.get("states", [])
    rows: list[dict[str, Any]] = []

    for state in states:
        if not isinstance(state, dict):
            continue

        yield_set, yield_fit, yield_err = normalize_param(state.get("yield"))
        centroid_set, centroid_fit, centroid_err = normalize_param(state.get("centroid_MeV"))
        width_set, width_fit, width_err = normalize_param(state.get("width_MeV"))

        rows.append(
            {
                "permutation": permutation,
                "fit": fit_number,
                "theta": theta,
                "state": state.get("index"),
                "shape": state.get("shape"),
                "yield_set": yield_set,
                "yield_fit": yield_fit,
                "yield_err": yield_err,
                "centroid_set": centroid_set,
                "centroid_fit": centroid_fit,
                "centroid_err": centroid_err,
                "width_set": width_set,
                "width_fit": width_fit,
                "width_err": width_err,
                "phase": phase,
                "phase_err": phase_err,
                "chi2": chi2,
                "ndf": ndf,
                "chi2_ndf": chi2_ndf,
            }
        )

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract BWFit state parameters from JSON summaries into a CSV."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="JSON files and/or directories containing JSON files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("state_fit_results.csv"),
        help="Output CSV path. Default: state_fit_results.csv",
    )
    parser.add_argument(
        "--permutation",
        default=None,
        help="Label for the permutation column. Default: output filename stem.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search directories recursively for JSON files.",
    )
    args = parser.parse_args()

    permutation = args.permutation or args.output.stem
    paths = json_paths(args.inputs, args.recursive)

    if not paths:
        raise SystemExit("No JSON files found.")

    all_rows: list[dict[str, Any]] = []
    for fit_number, path in enumerate(paths, start=1):
        rows = extract_rows(path, permutation, fit_number)
        if not rows:
            print(f"[warn] No states found in {path}")
            continue
        all_rows.extend(rows)
        print(f"[ok] {path}: extracted {len(rows)} states")

    if not all_rows:
        raise SystemExit("No state rows extracted.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\nWrote {len(all_rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
