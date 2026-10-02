#!/usr/bin/env python3
"""Plot BWFit state parameters versus spectrometer angle.

Examples
--------
python plot_fit_parameters_vs_angle.py fixed_phase_free_cent_width.csv
python plot_fit_parameters_vs_angle.py fixed_phase_free_cent_width.csv --show
python plot_fit_parameters_vs_angle.py fixed_phase_free_cent_width.csv -o plots --format pdf
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


REQUIRED_COLUMNS = {
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
    "chi2_ndf",
}


def load_results(csv_path: Path) -> pd.DataFrame:
    """Read and validate a BWFit parameter-summary CSV."""
    df = pd.read_csv(csv_path)

    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"CSV is missing required columns: {missing_text}")

    numeric_columns = [
        "theta",
        "state",
        "yield_set",
        "yield_fit",
        "yield_err",
        "centroid_set",
        "centroid_fit",
        "centroid_err",
        "width_set",
        "width_fit",
        "width_err",
        "chi2_ndf",
    ]
    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(subset=["theta", "state"]).copy()
    df["state"] = df["state"].astype(int)
    return df.sort_values(["state", "theta"])


def state_label(state_df: pd.DataFrame) -> str:
    """Create a compact legend label for one state."""
    state = int(state_df["state"].iloc[0])
    shape = str(state_df["shape"].iloc[0])
    return f"State {state} ({shape})"


def plot_parameter(
    df: pd.DataFrame,
    *,
    set_column: str,
    fit_column: str,
    error_column: str,
    ylabel: str,
    title: str,
    output_path: Path,
    show: bool,
    log_y: bool = False,
) -> None:
    """Plot set and fitted values for every state in one figure."""
    fig, ax = plt.subplots(figsize=(11, 8))

    for _, state_df in df.groupby("state", sort=True):
        state_df = state_df.sort_values("theta")
        label = state_label(state_df)

        # Fitted values with uncertainties.
        fit_container = ax.errorbar(
            state_df["theta"],
            state_df[fit_column],
            yerr=state_df[error_column],
            fmt="o-",
            capsize=4,
            linewidth=1.5,
            markersize=5,
            label=f"{label}: fit",
        )

        # Match the set-value curve to the corresponding fitted curve color.
        fit_color = fit_container.lines[0].get_color()
        ax.plot(
            state_df["theta"],
            state_df[set_column],
            "x--",
            linewidth=1.4,
            color=fit_color,
            label=f"{label}: set",
        )

    ax.set_xlabel(r"$\theta_{\mathrm{lab}}$ (deg)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)

    if log_y:
        ax.set_yscale("log")

    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_fit_quality(
    df: pd.DataFrame,
    *,
    output_path: Path,
    show: bool,
) -> None:
    """Plot one chi-square/NDF value per fit angle."""
    quality = (
        df[["theta", "chi2_ndf"]]
        .drop_duplicates(subset=["theta"])
        .sort_values("theta")
    )

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(quality["theta"], quality["chi2_ndf"], "o-")
    ax.set_xlabel(r"$\theta_{\mathrm{lab}}$ (deg)")
    ax.set_ylabel(r"$\chi^2 / \mathrm{NDF}$")
    ax.set_title("Fit quality versus angle")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot BWFit set and fitted state parameters versus angle."
    )
    parser.add_argument("csv", type=Path, help="Input fit-summary CSV file")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("fit_parameter_plots"),
        help="Directory for saved figures (default: fit_parameter_plots)",
    )
    parser.add_argument(
        "--format",
        choices=("png", "pdf", "svg"),
        default="png",
        help="Output figure format (default: png)",
    )
    parser.add_argument(
        "--log-yield",
        action="store_true",
        help="Use a logarithmic y-axis for the yield figure",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display figures interactively in addition to saving them",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = load_results(args.csv)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    suffix = args.format

    plot_parameter(
        df,
        set_column="centroid_set",
        fit_column="centroid_fit",
        error_column="centroid_err",
        ylabel="Centroid (MeV)",
        title="State centroids versus angle",
        output_path=args.output_dir / f"centroids_vs_angle.{suffix}",
        show=args.show,
    )

    plot_parameter(
        df,
        set_column="width_set",
        fit_column="width_fit",
        error_column="width_err",
        ylabel="Width (MeV)",
        title="State widths versus angle",
        output_path=args.output_dir / f"widths_vs_angle.{suffix}",
        show=args.show,
    )

    plot_parameter(
        df,
        set_column="yield_set",
        fit_column="yield_fit",
        error_column="yield_err",
        ylabel="Yield",
        title="State yields versus angle",
        output_path=args.output_dir / f"yields_vs_angle.{suffix}",
        show=args.show,
        log_y=args.log_yield,
    )

    plot_fit_quality(
        df,
        output_path=args.output_dir / f"chi2_ndf_vs_angle.{suffix}",
        show=args.show,
    )

    print(f"Saved plots to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError, pd.errors.ParserError) as exc:
        raise SystemExit(f"error: {exc}") from exc
