#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent

FIT_CONFIGS = {
    0: (
        "free_centroid_free_width.csv",
        "Free Centroid + Free Width"
    ),
    1: (
        "free_centroid_fixed_width.csv",
        "Free Centroid + Fixed Width"
    ),
    2: (
        "fixed_centroid_free_width.csv",
        "Fixed Centroid + Free Width"
    ),
    3: (
        "fixed_centroid_fixed_width.csv",
        "Fixed Centroid + Fixed Width"
    ),
}


def load_fit_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)

    # Replace this with real angles once known.
    # Example: df["theta"] = [10, 15, 20, 35, 40]
    if "theta" not in df.columns:
        df["theta"] = df["fit"]

    return df.sort_values("theta")


# def plot_2x2(df: pd.DataFrame, output: Path) -> None:
def plot_2x2(
    df: pd.DataFrame,
    output: Path,
    title: str
) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(13, 10), sharex=True)

    ax = axs[0, 0]
    ax.errorbar(df["theta"], df["trapped_centroid_fit"],
                yerr=df["trapped_centroid_err"],
                fmt="o-", capsize=4, label="Trapped fit")

    ax.plot(df["theta"], df["trapped_centroid_set"],
            "o--", label="Trapped set")

    ax.errorbar(df["theta"], df["superrad_centroid_fit"],
                yerr=df["superrad_centroid_err"],
                fmt="s-", capsize=4, label="Superrad fit")

    ax.plot(df["theta"], df["superrad_centroid_set"],
            "s--", label="Superrad set")

    ax.set_ylabel("Centroid (MeV)")
    ax.grid(alpha=0.3)
    # super extra legend handling because errorbar and plot don't play nice together
    handles, labels = ax.get_legend_handles_labels()
    legend_dict = dict(zip(labels, handles))
    ax.legend(
        [
            legend_dict["Trapped set"],
            legend_dict["Trapped fit"],
            legend_dict["Superrad set"],
            legend_dict["Superrad fit"],
        ],
        [
            "Trapped set",
            "Trapped fit",
            "Superrad set",
            "Superrad fit",
        ]
    )



    # (1,2) Width
    ax = axs[0, 1]
    ax.errorbar(df["theta"], df["trapped_width_fit"],
            yerr=df["trapped_width_err"],
            fmt="o-", capsize=4, label="Trapped fit")

    ax.plot(df["theta"], df["trapped_width_set"],
            "o--", label="Trapped set")
    
    ax.errorbar(df["theta"], df["superrad_width_fit"],
            yerr=df["superrad_width_err"],
            fmt="s-", capsize=4, label="Superrad fit")
    ax.plot(df["theta"], df["superrad_width_set"],
            "s--", label="Superrad set")    
    
    ax.set_ylabel("Width (MeV)")
    ax.grid(alpha=0.3)
        # super extra legend handling because errorbar and plot don't play nice together
    handles, labels = ax.get_legend_handles_labels()
    legend_dict = dict(zip(labels, handles))
    ax.legend(
        [
            legend_dict["Trapped set"],
            legend_dict["Trapped fit"],
            legend_dict["Superrad set"],
            legend_dict["Superrad fit"],
        ],
        [
            "Trapped set",
            "Trapped fit",
            "Superrad set",
            "Superrad fit",
        ]
    )
    


    # (2,1) Yield
    ax = axs[1, 0]
    ax.errorbar(df["theta"], df["trapped_yield_fit"],
                yerr=df["trapped_yield_err"],
                fmt="o-", capsize=4, label="Trapped fit")

    ax.plot(df["theta"], df["trapped_yield_set"],
        "o--", label="Trapped set")
    ax.errorbar(df["theta"], df["superrad_yield_fit"],
                yerr=df["superrad_yield_err"],
                fmt="s-", capsize=4, label="Superrad fit") 
    ax.plot(df["theta"], df["superrad_yield_set"],
        "s--", label="Superrad set")

    ax.set_yscale("log")
    ax.set_ylabel("Yield")
    ax.set_xlabel(r"$\theta_{\mathrm{lab}}$ (deg)")
    ax.grid(alpha=0.3)
        # super extra legend handling because errorbar and plot don't play nice together
    handles, labels = ax.get_legend_handles_labels()
    legend_dict = dict(zip(labels, handles))
    ax.legend(
        [
            legend_dict["Trapped set"],
            legend_dict["Trapped fit"],
            legend_dict["Superrad set"],
            legend_dict["Superrad fit"],
        ],
        [
            "Trapped set",
            "Trapped fit",
            "Superrad set",
            "Superrad fit",
        ]
    )

    # (2,2) Phase and chi2/ndf
    ax = axs[1, 1]
    ax.errorbar(df["theta"], df["phase"], yerr=df["phase_err"],
                fmt="o", capsize=4, label="Interference phase")
    # ax2 = ax.twinx()
    # ax2.plot(df["theta"], df["chi2_ndf"], "s--", label=r"$\chi^2/\mathrm{NDF}$")

    ax.set_ylabel("Phase (rad)")
    # ax2.set_ylabel(r"$\chi^2/\mathrm{NDF}$")
    ax.set_xlabel(r"$\theta_{\mathrm{lab}}$ (deg)")

    lines1, labels1 = ax.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    # ax.legend(lines1 + lines2, labels1 + labels2)
    ax.legend(lines1, labels1)
    ax.grid(alpha=0.3)

    fig.suptitle(
        f"{title} Breit-Wigner Fit Parameters",
        fontsize=16
    )
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.show()  # Show the figure
    plt.close(fig)



def main() -> None:
    fit_choice = 1

    csv_file, title = FIT_CONFIGS[fit_choice]

    csv_path = SCRIPT_DIR / csv_file
    output = SCRIPT_DIR / f"{Path(csv_file).stem}_2x2.png"

    print(f"Reading: {csv_path}")

    df = load_fit_csv(csv_path)
    plot_2x2(df, output, title)

    print(f"Saved: {output}")


if __name__ == "__main__":
    main()