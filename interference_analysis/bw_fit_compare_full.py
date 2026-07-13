#!/usr/bin/env python3

from pathlib import Path
import math
import pandas as pd
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent


def load_all_states_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    if "theta" not in df.columns:
        df["theta"] = df["fit"]

    return df.sort_values(["state", "theta"])


def plot_parameter_grid(
    df: pd.DataFrame,
    parameter: str,
    ylabel: str,
    output: Path,
    title: str,
    logy: bool = False,
) -> None:
    states = sorted(df["state"].unique())
    n_states = len(states)

    ncols = 2
    nrows = math.ceil(n_states / ncols)

    fig, axs = plt.subplots(
        nrows,
        ncols,
        figsize=(13, 3.2 * nrows),
        sharex=True,
    )

    axs = axs.flatten()

    set_col = f"{parameter}_set"
    fit_col = f"{parameter}_fit"
    err_col = f"{parameter}_err"

    for ax, state in zip(axs, states):
        sub = df[df["state"] == state]

        shape = sub["shape"].iloc[0]

        ax.plot(
            sub["theta"],
            sub[set_col],
            "o--",
            label="Set",
        )

        ax.errorbar(
            sub["theta"],
            sub[fit_col],
            yerr=sub[err_col],
            fmt="s-",
            capsize=4,
            label="Fit",
        )

        ax.set_title(f"State {state} ({shape})")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)

        if logy:
            ax.set_yscale("log")

        ax.legend()

    for ax in axs[n_states:]:
        ax.axis("off")

    for ax in axs[-ncols:]:
        ax.set_xlabel(r"$\theta_{\mathrm{lab}}$ (deg)")

    fig.suptitle(title, fontsize=16)
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.show()
    plt.close(fig)


def main() -> None:
    csv_path = SCRIPT_DIR / "free_centroid_free_width_all_states.csv"
    df = load_all_states_csv(csv_path)

    plot_parameter_grid(
        df,
        parameter="centroid",
        ylabel="Centroid (MeV)",
        output=SCRIPT_DIR / "all_states_centroids.png",
        title="All-State Centroids: Set vs Fit",
    )

    plot_parameter_grid(
        df,
        parameter="width",
        ylabel="Width (MeV)",
        output=SCRIPT_DIR / "all_states_widths.png",
        title="All-State Widths: Set vs Fit",
    )

    plot_parameter_grid(
        df,
        parameter="yield",
        ylabel="Yield",
        output=SCRIPT_DIR / "all_states_yields.png",
        title="All-State Yields: Set vs Fit",
        logy=True,
    )


if __name__ == "__main__":
    main()