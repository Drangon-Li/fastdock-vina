#!/usr/bin/env python3
"""Render the README performance figure from the published benchmark CSV.

Requires Python 3.8+ and Matplotlib (python3 -m pip install matplotlib).
Run from any directory:
    python3 tools/plot_speedup.py
"""
from pathlib import Path
import csv
import json
import math
import statistics
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "fastdock-matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

ROOT = Path(__file__).resolve().parents[1]
TEAL = "#087f8c"
ORANGE = "#cf731d"
INK = "#26343b"
MUTED = "#5d6a70"


def load_data():
    with (ROOT / "benchmarks/benchmark.csv").open(newline="") as source:
        rows = list(csv.DictReader(source))
    assert len(rows) == 25
    for row in rows:
        assert int(row["repeats_per_variant"]) == 3
        for variant in ("baseline", "v2"):
            samples = json.loads(row[f"{variant}_wall_seconds_raw"])
            assert len(samples) == 3
            median = statistics.median(samples)
            assert math.isclose(median, float(row[f"{variant}_wall_seconds"]), rel_tol=1e-12)
            row[variant] = median
        row["speedup"] = row["baseline"] / row["v2"]
    return rows


def case_label(case):
    labels = {
        "basic-vina-dock-cpu1": "Basic / Vina / 1 thread",
        "basic-vina-dock-seed20261013": "Basic / Vina / alternate seed",
    }
    if case in labels:
        return labels[case]
    parts = case.split("-")
    return " / ".join(p.upper() if p == "ad4" else p.capitalize() for p in parts if p != "dock")


def panel(ax, rows, color, title, xmax, ticks):
    y = list(range(len(rows)))
    values = [r["speedup"] for r in rows]
    ax.axvline(1, color="#839197", linewidth=1.3, linestyle=(0, (3, 3)), zorder=1)
    ax.hlines(y, 1, values, color=color, linewidth=3, alpha=0.75, zorder=2)
    ax.scatter(values, y, s=49, color=color, edgecolors="white", linewidths=0.8, zorder=3)
    for i, value in enumerate(values):
        ax.annotate(f"{value:.3f}×", (value, i), xytext=(8, 0), textcoords="offset points",
                    va="center", fontsize=11, color=INK, weight="normal")
    ax.set_yticks(y, [case_label(r["case"]) for r in rows], fontsize=11, color=INK)
    ax.set_ylim(len(rows) - 0.45, -0.65)
    ax.set_xlim(0.99 if xmax < 1.5 else 0.95, xmax)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}×"))
    ax.tick_params(axis="both", length=0, pad=8, labelcolor=MUTED)
    ax.grid(axis="x", color="#e7edef", linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, loc="left", fontsize=15, color=INK, weight="bold", pad=19)
    ax.set_xlabel("Wall-time speedup · baseline / FastDock Vina", color=MUTED, fontsize=11, labelpad=14)


def main():
    rows = load_data()
    full = [r for r in rows if r["group"] == "full-search"]
    local = [r for r in rows if r["group"] == "score/local"]
    assert (len(full), len(local)) == (19, 6)
    geometric = statistics.geometric_mean(r["speedup"] for r in rows)
    baseline_total = sum(r["baseline"] for r in full)
    fastdock_total = sum(r["v2"] for r in full)
    reduction = (1 - fastdock_total / baseline_total) * 100
    aggregate = baseline_total / fastdock_total
    assert math.isclose(geometric, 1.2283158586567096, rel_tol=1e-12)
    assert math.isclose(reduction, 7.811626341647571, rel_tol=1e-12)
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "path", "font.size": 12})
    fig = plt.figure(figsize=(14, 11.2), facecolor="white")
    fig.text(0.055, 0.955, "FastDock Vina performance", fontsize=25, color=INK, weight="bold")
    fig.text(0.055, 0.922, "CPU-only  ·  25 configurations  ·  median wall time of n = 3 runs per variant", fontsize=12.5, color=MUTED)
    fig.text(0.055, 0.868, f"{geometric:.3f}×", fontsize=28, color=TEAL, weight="bold")
    fig.text(0.055, 0.839, "Equal-weight geometric mean", fontsize=13, color=INK, weight="bold")
    fig.text(0.055, 0.815, "All 25 configurations, including score/local", fontsize=11.5, color=MUTED)
    fig.text(0.55, 0.868, f"{reduction:.2f}% less time", fontsize=28, color=TEAL, weight="bold")
    fig.text(0.55, 0.839, "Full docking: summed median wall times", fontsize=13, color=INK, weight="bold")
    fig.text(0.55, 0.815, f"19 configurations  ·  {aggregate:.3f}× aggregate speedup", fontsize=11.5, color=MUTED)
    fig.add_artist(plt.Line2D([0.055, 0.955], [0.79, 0.79], transform=fig.transFigure, color="#dce5e8", linewidth=1))
    left = fig.add_axes([0.25, 0.105, 0.255, 0.625])
    right = fig.add_axes([0.735, 0.475, 0.215, 0.255])
    panel(left, full, TEAL, "Full docking · 19 configurations", 1.285, [1, 1.1, 1.2])
    panel(right, local, ORANGE, "Score / local · 6 configurations", 2.58, [1, 1.5, 2])
    fig.text(0.61, 0.345, "Reading the figure", color=INK, fontsize=13, weight="bold")
    fig.text(0.61, 0.311, "1× = same wall time as baseline\nGreater than 1× = faster\n\nPanels use different horizontal scales.\nScore/local runs are specialized operations;\nfull docking is reported separately above.", color=MUTED, fontsize=11.5, linespacing=1.7, va="top")
    fig.text(0.055, 0.035, "Source: benchmarks/benchmark.csv  ·  Ratios computed from per-configuration medians  ·  No error bars shown", color=MUTED, fontsize=10.5)
    output = ROOT / "docs/figures"
    output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(output / f"speedup-summary.{extension}", dpi=160, facecolor="white", metadata={"Creator": "FastDock Vina benchmark plotting script"})
    plt.close(fig)
    print(f"25-config geometric mean: {geometric:.16f}x")
    print(f"19 full-docking totals: baseline={baseline_total:.6f}s, FastDock Vina={fastdock_total:.6f}s")
    print(f"19 full-docking aggregate: {aggregate:.16f}x; wall-time reduction={reduction:.15f}%")
    print("Wrote docs/figures/speedup-summary.png and docs/figures/speedup-summary.svg")


if __name__ == "__main__":
    main()
