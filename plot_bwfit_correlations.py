#!/usr/bin/env python3
"""
Plot BWFit correlation diagnostics from a JSON fit summary.

Creates:
  1. Free-parameter correlation heatmap
  2. Correlations with interference phase
  3. Local covariance confidence contours for trapped centroid vs phase

Important:
The third plot is a LOCAL Gaussian approximation derived from the covariance
matrix. It is not a profiled 2D chi-square scan.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("json_file")
    parser.add_argument("--outdir", default="correlation_plots")
    args = parser.parse_args()

    with open(args.json_file, "r") as f:
        data = json.load(f)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fit = data["fit_result"]
    labels = fit["parameter_labels"]
    fixed = fit["parameter_fixed"]
    corr = np.array(fit["correlation_matrix"], dtype=float)
    cov = np.array(fit["covariance_matrix"], dtype=float)

    pretty = {
        "State0_Yield": r"$Y_{6.875}$",
        "State1_Yield": r"$Y_{7.4665}$",
        "State2_Yield": r"$Y_{7.524}$",
        "State3_Yield": r"$Y_{\mathrm{T}}$",
        "State3_M": r"$E_{\mathrm{T}}$",
        "State3_G": r"$\Gamma_{\mathrm{T}}$",
        "State4_Yield": r"$Y_{\mathrm{SR}}$",
        "State4_M": r"$E_{\mathrm{SR}}$",
        "State4_G": r"$\Gamma_{\mathrm{SR}}$",
        "State5_Yield": r"$Y_{9.498}$",
        "State5_M": r"$E_{9.498}$",
        "State5_G": r"$\Gamma_{9.498}$",
        "Phase_delta": r"$\delta$",
        "BG_A0": r"$A_0$",
        "BG_A1": r"$A_1$",
        "BG_E0": r"$E_0$",
    }
    plot_labels = [pretty.get(x, x) for x in labels]
    free_idx = [i for i, is_fixed in enumerate(fixed) if not is_fixed]
    free_labels = [plot_labels[i] for i in free_idx]
    corr_free = corr[np.ix_(free_idx, free_idx)]

    # Heatmap
    fig, ax = plt.subplots(figsize=(11, 9))
    im = ax.imshow(corr_free, vmin=-1, vmax=1, aspect="equal")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"Correlation coefficient $\rho$")
    ax.set_xticks(range(len(free_labels)))
    ax.set_yticks(range(len(free_labels)))
    ax.set_xticklabels(free_labels, rotation=45, ha="right")
    ax.set_yticklabels(free_labels)
    ax.set_title(r"Correlation matrix (free parameters only)")
    for i in range(len(free_labels)):
        for j in range(len(free_labels)):
            ax.text(j, i, f"{corr_free[i,j]:.2f}",
                    ha="center", va="center", fontsize=7)
    fig.tight_layout()
    fig.savefig(outdir / "correlation_heatmap.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "correlation_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)




    # Phase correlations
    phase_idx = labels.index("Phase_delta")
    phase_pairs = []
    for i in free_idx:
        if i == phase_idx:
            continue
        phase_pairs.append((plot_labels[i], corr[phase_idx, i]))
    phase_pairs.sort(key=lambda t: abs(t[1]), reverse=True)

    names = [x[0] for x in phase_pairs][::-1]
    vals = [x[1] for x in phase_pairs][::-1]

    fig, ax = plt.subplots(figsize=(8, 7))
    y = np.arange(len(names))
    ax.barh(y, vals)
    ax.axvline(0, linewidth=1)
    ax.axvline(0.7, linestyle="--", linewidth=1)
    ax.axvline(-0.7, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_xlim(-1, 1)
    ax.set_xlabel(r"Correlation with phase, $\rho(\delta,p)$")
    ax.set_title(r"Correlations with interference phase")
    fig.tight_layout()
    fig.savefig(outdir / "phase_correlations.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "phase_correlations.pdf", bbox_inches="tight")
    plt.close(fig)




    # Local covariance contours: trapped centroid vs phase
    et_idx = labels.index("State3_M")
    x0 = data["fit_result"]["states"][3]["centroid_MeV"]["fit"]
    y0 = data["fit_result"]["interference"]["phase_fit_rad"]

    subcov = cov[np.ix_([et_idx, phase_idx], [et_idx, phase_idx])]
    eigvals, eigvecs = np.linalg.eigh(subcov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    angle_deg = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))

    # levels = [(2.3, r"$\Delta \chi^2 = 2.30$"), (6.18, r"$\Delta \chi^2 = 6.18$"), (11.83, r"$\Delta \chi^2 = 11.83$")]
    levels = [
    (2.30, r"$\Delta \chi^2 = 2.30$", "-"),
    (6.18, r"$\Delta \chi^2 = 6.18$", "--"),
    (11.83, r"$\Delta \chi^2 = 11.83$", ":"),
]

    fig, ax = plt.subplots(figsize=(8, 6))
    # for dchi2, label in levels:
    for dchi2, label, linestyle in levels:
        width = 2 * math.sqrt(dchi2 * eigvals[0])
        height = 2 * math.sqrt(dchi2 * eigvals[1])
        ell = Ellipse(
            (x0, y0),
            width=width,
            height=height,
            angle=angle_deg,
            fill=False,
            linewidth=1.8,
            linestyle=linestyle,
            label=f"{label} local covariance contour",
        )
        ax.add_patch(ell)

    ax.scatter([x0], [y0], s=35, zorder=3)
    ax.axvline(x0, linestyle=":", linewidth=1)
    ax.axhline(y0, linestyle=":", linewidth=1)
    ax.set_xlabel(r"Trapped-state centroid $E_{\mathrm{T}}$ (MeV)", fontsize=12)
    ax.set_ylabel(r"Interference phase $\delta$ (rad)", fontsize=12)
    ax.set_title(
        # r"Local $E_{\mathrm{T}}$--$\delta$ covariance contours"
        # + "\n +"
        rf"$\rho = {corr[et_idx, phase_idx]:.3f}$"
    )
    ax.legend()
    ax.autoscale()
    fig.tight_layout()
    fig.savefig(outdir / "phase_vs_trapped_centroid_covariance.png",
                dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "phase_vs_trapped_centroid_covariance.pdf",
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
