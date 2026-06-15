"""
Plot COCOA.jl computational performance scaling.

Produces separate linear-scale figures:
    runtime_scaling_linear (runtime in hours)
    memory_scaling_linear (memory usage in GB)

Color-coded by variant (enzyme binding mechanism), with per-variant
Theil-Sen median lines + shaded IQR bands (25th--75th percentile).
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats as sp_stats
import statsmodels.regression.quantile_regression as smq
import statsmodels.api as sm
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import style as thesis_style
thesis_style.apply()

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_PATH = REPO_ROOT / "results" / "kinetic" / "resource_analysis.csv"
OUT_DIR   = REPO_ROOT / "figures"
OUT_DIR.mkdir(exist_ok=True, parents=True)

VARIANT_ORDER = [
    "no_split",
    "random_0",
    "random_25",
    "random_50",
    "random_75",
    "random_100",
]

VARIANT_LABELS = {
    "no_split": "No splitting",
    "random_0": "0%",
    "random_25": "25%",
    "random_50": "50%",
    "random_75": "75%",
    "random_100": "100%",
}


def fit_theilsen(log_x, log_y):
    """Theil-Sen median slope in log-log space. Returns (slope, intercept)."""
    result = sp_stats.theilslopes(log_y, log_x)
    return float(result.slope), float(result.intercept)


def fit_quantile(log_x, log_y, q):
    """Quantile regression at quantile q in log-log space. Returns (slope, intercept)."""
    X = sm.add_constant(log_x)
    model = smq.QuantReg(log_y, X)
    res = model.fit(q=q, max_iter=2000)
    intercept, slope = float(res.params[0]), float(res.params[1])
    return slope, intercept


def predict_line(slope, intercept, x_range):
    """Convert log-log fit back to original scale."""
    return 10 ** (intercept + slope * np.log10(x_range))


def make_single_figure(df_ok, y_col, y_label, variant_colors, xscale="linear", yscale="linear"):
    FS_LABEL  = thesis_style.FS_LABEL
    FS_TICK   = thesis_style.FS_TICK
    FS_LEGEND = thesis_style.FS_LEGEND
    FS_LTITLE = thesis_style.FS_LTITLE
    fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))

    x_min = df_ok["n_complexes"].min()
    x_max = df_ok["n_complexes"].max()

    for variant in VARIANT_ORDER:
        sub = df_ok[df_ok["variant"] == variant]
        if variant not in df_ok["variant"].values or len(sub) < 5:
            continue

        color = variant_colors[variant]
        log_x = np.log10(sub["n_complexes"].values)
        log_y = np.log10(sub[y_col].values)

        ax.scatter(
            sub["n_complexes"],
            sub[y_col],
            alpha=0.40,
            s=25,
            color=color,
            label=VARIANT_LABELS[variant],
            zorder=2,
        )

        vx_min = sub["n_complexes"].min()
        vx_max = sub["n_complexes"].max()
        vx_range = np.logspace(np.log10(vx_min), np.log10(vx_max), 200)

        slope_med, int_med = fit_theilsen(log_x, log_y)
        ax.plot(
            vx_range,
            predict_line(slope_med, int_med, vx_range),
            color=color,
            linewidth=2.5,
            zorder=3,
        )

        slope_lo, int_lo = fit_quantile(log_x, log_y, 0.25)
        slope_hi, int_hi = fit_quantile(log_x, log_y, 0.75)
        ax.fill_between(
            vx_range,
            predict_line(slope_lo, int_lo, vx_range),
            predict_line(slope_hi, int_hi, vx_range),
            color=color,
            alpha=0.18,
            zorder=1,
        )

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    if xscale == "log":
        # Explicitly place labels at sensible values across the data range
        x_ticks = [t for t in [3000, 4000, 5000, 6000, 8000, 10000, 15000, 20000, 30000, 40000, 50000]
                   if x_min * 0.9 <= t <= x_max * 1.1]
        ax.xaxis.set_major_locator(ticker.FixedLocator(x_ticks))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    else:
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    ax.set_xlabel("Number of Complexes", fontsize=FS_LABEL)
    ax.set_ylabel(y_label, fontsize=FS_LABEL)
    ax.tick_params(axis="both", labelsize=FS_TICK)
    ax.tick_params(axis="x", labelrotation=30)
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(
        title="Random binding",
        fontsize=FS_LEGEND,
        title_fontsize=FS_LTITLE,
        framealpha=0.85,
        edgecolor="0.7",
    )

    plt.tight_layout()
    return fig


def main():
    df = pd.read_csv(DATA_PATH)
    df_ok = df[df["job_status"] == "SUCCESS"].copy()
    df_ok["duration_hours"] = df_ok["duration_minutes"] / 60.0

    if "slurm_max_rss_gb" in df_ok.columns and df_ok["slurm_max_rss_gb"].notna().any():
        memory_col = "slurm_max_rss_gb"
    else:
        memory_col = "memory_allocated_gb"

    colors = [thesis_style.PALETTE[i] for i in range(len(VARIANT_ORDER))]
    variant_colors = dict(zip(VARIANT_ORDER, colors))

    outputs = [
        ("duration_hours", "Runtime (hours)", "runtime_scaling_linear"),
        (memory_col, "Memory Usage (GB)", "memory_scaling_linear"),
    ]

    for y_col, y_label, stem in outputs:
        fig = make_single_figure(
            df_ok,
            y_col,
            y_label,
            variant_colors,
            xscale="linear",
            yscale="linear",
        )
        fig.savefig(OUT_DIR / f"{stem}.pdf", dpi=300, bbox_inches="tight")
        fig.savefig(OUT_DIR / f"{stem}.svg", bbox_inches="tight")
        fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
        plt.close()

    print(f"Figures saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
