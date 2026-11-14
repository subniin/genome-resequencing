"""Runs experiment steps 0-4: generates data, assembles, evaluates, and records time and memory.

Usage: python scripts/run_experiments.py [--steps 0 1 2 3 4] [--skip-greedy-100k]
Writes results/step{N}.csv. Build artifacts go to build/, data to data/.
"""
import argparse
import csv
import math
import os
import subprocess
import sys
import time

import psutil

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BUILD = os.path.join(ROOT, "build")
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
SEED = 1
D = 2
K = 14

TOOLS = ["ssake", "greedy_overlap"]
COLUMNS = [
    "step", "variable", "value", "tool", "N", "L", "coverage", "snps", "o", "M",
    "seconds", "peak_mem_mb", "accuracy", "snp_recall", "snp_precision",
    "n_contigs", "n50", "covered_frac", "misassemblies", "unaligned_contigs",
]


def exe(name):
    return os.path.join(BUILD, name + ".exe")


def build_tools():
    os.makedirs(BUILD, exist_ok=True)
    for src in ["gen_reference", "gen_reads", "ssake", "greedy_overlap", "evaluate"]:
        subprocess.run(
            ["g++", "-O2", "-std=c++17", "-static",
             os.path.join(ROOT, "src", src + ".cpp"), "-o", exe(src)],
            check=True,
        )


def run_measured(args):
    """Runs a command; returns (wall seconds, peak RSS of the process tree in MB, stdout)."""
    t0 = time.perf_counter()
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    ps = psutil.Process(proc.pid)
    peak = 0
    while proc.poll() is None:
        try:
            procs = [ps] + ps.children(recursive=True)
            rss = sum(p.memory_info().rss for p in procs)
            peak = max(peak, rss)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        time.sleep(0.02)
    out, err = proc.communicate()
    seconds = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"{args[0]} failed: {err.strip()}")
    return seconds, peak / (1024 * 1024), out


def ensure_reference(N):
    path = os.path.join(DATA, f"ref_{N}.txt")
    if not os.path.exists(path):
        os.makedirs(DATA, exist_ok=True)
        subprocess.run([exe("gen_reference"), "--N", str(N), "--k", str(K),
                        "--seed", str(SEED), "--out", path], check=True)
    return path


def run_case(step, variable, value, tool, N, L, coverage, snp_count, o):
    """One assembly + evaluation. Returns a CSV row dict."""
    ref = ensure_reference(N)
    tag = f"N{N}_L{L}_c{coverage}_s{snp_count}_o{o}"
    run_dir = os.path.join(DATA, "runs", tag)
    os.makedirs(run_dir, exist_ok=True)

    subprocess.run([exe("gen_reads"), "--ref", ref, "--L", str(L), "--coverage", str(coverage),
                    "--snps", str(snp_count), "--seed", str(SEED), "--out", run_dir], check=True)
    reads = os.path.join(run_dir, "reads.txt")
    contigs = os.path.join(run_dir, f"contigs_{tool}.fa")
    seconds, peak_mb, _ = run_measured(
        [exe(tool), "--reads", reads, "--o", str(o), "--D", str(D), "--out", contigs])

    out = subprocess.run(
        [exe("evaluate"), "--contigs", contigs, "--ref", ref,
         "--sample", os.path.join(run_dir, "sample.txt"),
         "--truth", os.path.join(run_dir, "truth_snps.tsv"), "--k", str(K)],
        check=True, capture_output=True, text=True).stdout
    metrics = dict(line.split("=", 1) for line in out.strip().splitlines())

    M = sum(1 for _ in open(reads))
    row = {
        "step": step, "variable": variable, "value": value, "tool": tool,
        "N": N, "L": L, "coverage": coverage, "snps": snp_count, "o": o, "M": M,
        "seconds": round(seconds, 3), "peak_mem_mb": round(peak_mb, 1),
    }
    row.update(metrics)
    print(f"[step {step}] {tool} {tag}: {seconds:.2f}s {peak_mb:.0f}MB "
          f"contigs={metrics['n_contigs']} recall={metrics['snp_recall']}")
    return row


def write_csv(step, rows):
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, f"step{step}.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("wrote", path)


def default_o(L):
    return math.ceil(L / 2)


def step0():
    return [run_case(0, "default", "-", t, 10000, 30, 10, 50, default_o(30)) for t in TOOLS]


def step1(skip_greedy_100k):
    rows = []
    for N in [10000, 100000, 1000000, 10000000]:
        rows.append(run_case(1, "N", N, "ssake", N, 30, 10, round(N * 0.005), default_o(30)))
        if N <= 10000 or (N == 100000 and not skip_greedy_100k):
            rows.append(run_case(1, "N", N, "greedy_overlap", N, 30, 10, round(N * 0.005), default_o(30)))
    return rows


def step2():
    rows = []
    for cov in [5, 10, 15]:
        for t in TOOLS:
            rows.append(run_case(2, "coverage", cov, t, 10000, 30, cov, 50, default_o(30)))
    return rows


def step3():
    rows = []
    for ratio, count in [(0.001, 10), (0.005, 50), (0.01, 100)]:
        for t in TOOLS:
            rows.append(run_case(3, "snp_ratio", ratio, t, 10000, 30, 10, count, default_o(30)))
    return rows


def step4():
    rows = []
    for L in [15, 30, 45]:
        for t in TOOLS:
            rows.append(run_case(4, "L", L, t, 10000, L, 10, 50, default_o(L)))
    for o in [10, 15, 20]:
        for t in TOOLS:
            rows.append(run_case(4, "o", o, t, 10000, 30, 10, 50, o))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    parser.add_argument("--skip-greedy-100k", action="store_true")
    args = parser.parse_args()

    build_tools()
    steps = {0: step0, 2: step2, 3: step3, 4: step4}
    for s in args.steps:
        if s == 1:
            rows = step1(args.skip_greedy_100k)
        else:
            rows = steps[s]()
        write_csv(s, rows)


if __name__ == "__main__":
    sys.exit(main())
