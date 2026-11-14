"""Plots results/step*.csv into results/*.png."""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def load(step):
    with open(os.path.join(RESULTS, f"step{step}.csv"), newline="") as f:
        return list(csv.DictReader(f))


def series(rows, tool, xkey, ykey):
    pts = sorted((float(r[xkey]), float(r[ykey])) for r in rows if r["tool"] == tool)
    return [p[0] for p in pts], [p[1] for p in pts]


def lander_waterman(cov, N=10000, L=30, o=15):
    lam = cov / L
    return N * lam * math.exp(-lam * (L - o))


def plot_step1():
    rows = load(1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for tool in ["ssake", "greedy_overlap"]:
        x, y = series(rows, tool, "N", "seconds")
        axes[0].plot(x, y, marker="o", label=tool)
        x, y = series(rows, tool, "N", "peak_mem_mb")
        axes[1].plot(x, y, marker="o", label=tool)
    for ax, ylabel in zip(axes, ["seconds", "peak memory (MB)"]):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("N (reference length)")
        ax.set_ylabel(ylabel)
        ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "step1_time_memory.png"), dpi=120)


def plot_step2():
    rows = load(2)
    fig, ax = plt.subplots(figsize=(6, 4))
    for tool in ["ssake", "greedy_overlap"]:
        x, y = series(rows, tool, "coverage", "n_contigs")
        ax.plot(x, y, marker="o", label=tool)
    cov = [5, 10, 15]
    ax.plot(cov, [lander_waterman(c) for c in cov], "k--", label="Lander-Waterman estimate")
    ax.set_yscale("log")
    ax.set_xlabel("coverage")
    ax.set_ylabel("contigs")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "step2_contigs_vs_coverage.png"), dpi=120)


def plot_step3():
    rows = load(3)
    fig, ax = plt.subplots(figsize=(6, 4))
    for tool in ["ssake", "greedy_overlap"]:
        x, y = series(rows, tool, "value", "snp_recall")
        ax.plot([v * 100 for v in x], y, marker="o", label=tool)
    ax.set_xlabel("SNP ratio (%)")
    ax.set_ylabel("SNP recall")
    ax.set_ylim(0, 1.05)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "step3_recall_vs_snp.png"), dpi=120)


if __name__ == "__main__":
    for s in [1, 2, 3]:
        if not os.path.exists(os.path.join(RESULTS, f"step{s}.csv")):
            print(f"skip step{s}: no csv")
    if os.path.exists(os.path.join(RESULTS, "step1.csv")):
        plot_step1()
    if os.path.exists(os.path.join(RESULTS, "step2.csv")):
        plot_step2()
    if os.path.exists(os.path.join(RESULTS, "step3.csv")):
        plot_step3()
    print("plots written to", os.path.abspath(RESULTS))
