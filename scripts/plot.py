"""Plots results/step*.csv into results/*.png."""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
TOOLS = ["ssake", "greedy_overlap"]


def load(step):
    with open(os.path.join(RESULTS, f"step{step}.csv"), newline="") as f:
        return list(csv.DictReader(f))


def series(rows, tool, xkey, ykey, variable=None):
    pts = sorted((float(r[xkey]), float(r[ykey])) for r in rows
                 if r["tool"] == tool and (variable is None or r["variable"] == variable))
    return [p[0] for p in pts], [p[1] for p in pts]


def lander_waterman(cov, N=100000, L=50, o=25):
    lam = cov / L
    return N * lam * math.exp(-lam * (L - o))


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, name), dpi=120)
    plt.close(fig)


def plot_step1():
    rows = load(1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for tool in TOOLS:
        x, y = series(rows, tool, "M", "seconds")
        axes[0].plot(x, y, marker="o", label=tool)
        x, y = series(rows, tool, "M", "peak_mem_mb")
        axes[1].plot(x, y, marker="o", label=tool)
    axes[1].axhline(2000, color="k", linestyle="--", label="target 2000 MB at M=2e7")
    for ax, ylabel in zip(axes, ["seconds", "peak memory (MB)"]):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("M (number of reads, L=50)")
        ax.set_ylabel(ylabel)
        ax.legend()
    save(fig, "step1_time_memory.png")


def plot_step2():
    rows = load(2)
    fig, ax = plt.subplots(figsize=(6, 4))
    for tool in TOOLS:
        x, y = series(rows, tool, "coverage", "n_contigs")
        ax.plot(x, y, marker="o", label=tool)
    cov = [5, 10, 15]
    ax.plot(cov, [lander_waterman(c) for c in cov], "k--", label="Lander-Waterman estimate")
    ax.set_yscale("log")
    ax.set_xlabel("coverage")
    ax.set_ylabel("contigs")
    ax.legend()
    save(fig, "step2_contigs_vs_coverage.png")


def plot_step3():
    rows = load(3)
    fig, ax = plt.subplots(figsize=(6, 4))
    for tool in TOOLS:
        x, y = series(rows, tool, "value", "snp_recall")
        ax.plot([v * 100 for v in x], y, marker="o", label=tool)
    ax.set_xlabel("SNP ratio (%)")
    ax.set_ylabel("SNP recall")
    ax.set_ylim(0, 1.05)
    ax.legend()
    save(fig, "step3_recall_vs_snp.png")


def plot_step4():
    rows = load(4)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for tool in TOOLS:
        x, y = series(rows, tool, "value", "n_contigs", variable="L")
        axes[0].plot(x, y, marker="o", label=tool)
        x, y = series(rows, tool, "value", "n50", variable="L")
        axes[1].plot(x, y, marker="o", label=tool)
    for ax, ylabel in zip(axes, ["contigs", "N50 (bp)"]):
        ax.set_xlabel("read length L (o = ceil(L/2))")
        ax.set_ylabel(ylabel)
        ax.set_yscale("log")
        ax.legend()
    save(fig, "step4_L_contigs_n50.png")


def plot_step5():
    rows = load(5)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for tool in TOOLS:
        x, y = series(rows, tool, "value", "accuracy", variable="err")
        axes[0].plot(x, y, marker="o", label=f"{tool} accuracy")
        x, y = series(rows, tool, "value", "snp_precision", variable="err")
        axes[0].plot(x, y, marker="s", linestyle=":", label=f"{tool} SNP precision")
        x, y = series(rows, tool, "value", "snp_recall", variable="D")
        axes[1].plot(x, y, marker="o", label=f"{tool} SNP recall")
    x, y = series(rows, "ssake_vote", "value", "accuracy", variable="err")
    axes[0].plot(x, y, marker="o", linestyle="--", label="ssake_vote accuracy")
    x, y = series(rows, "ssake_vote", "value", "snp_precision", variable="err")
    axes[0].plot(x, y, marker="s", linestyle="-.", label="ssake_vote SNP precision")
    axes[0].set_xlabel("sequencing error rate (D=2)")
    axes[0].set_ylabel("fraction")
    axes[1].set_xlabel("D (mismatches allowed, err=0.01)")
    axes[1].set_ylabel("SNP recall")
    for ax in axes:
        ax.legend(fontsize=8)
    save(fig, "step5_err_accuracy.png")


def plot_step6():
    rows = load(6)
    fig, ax = plt.subplots(figsize=(6, 4))
    for tool in TOOLS:
        x, y = series(rows, tool, "value", "misassemblies")
        ax.plot(x, y, marker="o", label=tool)
    ax.set_xlabel("repeat length (20 copies, L=50)")
    ax.set_ylabel("misassemblies")
    ax.legend()
    save(fig, "step6_repeat_misassembly.png")


if __name__ == "__main__":
    plotters = {1: plot_step1, 2: plot_step2, 3: plot_step3,
                4: plot_step4, 5: plot_step5, 6: plot_step6}
    for s, fn in plotters.items():
        if os.path.exists(os.path.join(RESULTS, f"step{s}.csv")):
            fn()
        else:
            print(f"skip step{s}: no csv")
    print("plots written to", os.path.abspath(RESULTS))
