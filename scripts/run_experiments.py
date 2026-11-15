"""Runs experiment steps 0-6 for the design-project conditions: generates data, assembles,
evaluates, and records time and memory.

Usage: python -u scripts/run_experiments.py [--steps 0 1 2 3 4 5 6]
Writes results/step{N}.csv. Build artifacts go to build/, data to data/.

Base case (step 0): N=100000, L=50, coverage 10, SNP 0.1%, o=ceil(L/2), D=2.
Step 1 keeps the base case but sets N = M*L/coverage, so M = 2e4 .. 2e7 (N = 1e5 .. 1e8).
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
K_GEN = 16    # k-mer length that makes the reference unique
K_EVAL = 24   # anchor k-mer length in the evaluator
L_DEFAULT = 50
COV_DEFAULT = 10
SNP_RATIO = 0.001
N_BASE = 100000
M_STEP1 = [20000, 200000, 2000000, 20000000]

TOOLS = ["ssake", "greedy_overlap"]
COLUMNS = [
    "step", "variable", "value", "tool", "N", "L", "coverage", "snps", "o", "D", "err",
    "repeats", "repeat_len", "M",
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


def ensure_reference(N, repeats, repeat_len):
    path = os.path.join(DATA, f"ref_N{N}_r{repeats}x{repeat_len}.txt")
    if not os.path.exists(path):
        os.makedirs(DATA, exist_ok=True)
        subprocess.run([exe("gen_reference"), "--N", str(N), "--k", str(K_GEN),
                        "--seed", str(SEED), "--repeats", str(repeats),
                        "--repeat-len", str(repeat_len), "--out", path], check=True)
    return path


def run_case(step, variable, value, tool, N, L, coverage, snp_count, o,
             D=D, err=0.0, repeats=0, repeat_len=0):
    """One assembly + evaluation. Returns a CSV row dict."""
    ref = ensure_reference(N, repeats, repeat_len)
    tag = f"N{N}_L{L}_c{coverage}_s{snp_count}_o{o}_D{D}_e{err}_r{repeats}x{repeat_len}"
    run_dir = os.path.join(DATA, "runs", tag)
    os.makedirs(run_dir, exist_ok=True)

    subprocess.run([exe("gen_reads"), "--ref", ref, "--L", str(L), "--coverage", str(coverage),
                    "--snps", str(snp_count), "--err", str(err), "--seed", str(SEED),
                    "--out", run_dir], check=True)
    reads = os.path.join(run_dir, "reads.txt")
    contigs = os.path.join(run_dir, f"contigs_{tool}.fa")
    seconds, peak_mb, _ = run_measured(
        [exe(tool), "--reads", reads, "--o", str(o), "--D", str(D), "--out", contigs])

    out = subprocess.run(
        [exe("evaluate"), "--contigs", contigs, "--ref", ref,
         "--sample", os.path.join(run_dir, "sample.txt"),
         "--truth", os.path.join(run_dir, "truth_snps.tsv"), "--k", str(K_EVAL)],
        check=True, capture_output=True, text=True).stdout
    metrics = dict(line.split("=", 1) for line in out.strip().splitlines())

    M = sum(1 for _ in open(reads))
    row = {
        "step": step, "variable": variable, "value": value, "tool": tool,
        "N": N, "L": L, "coverage": coverage, "snps": snp_count, "o": o, "D": D,
        "err": err, "repeats": repeats, "repeat_len": repeat_len, "M": M,
        "seconds": round(seconds, 3), "peak_mem_mb": round(peak_mb, 1),
    }
    row.update(metrics)
    print(f"[step {step}] {tool} {tag}: {seconds:.2f}s {peak_mb:.0f}MB "
          f"contigs={metrics['n_contigs']} recall={metrics['snp_recall']}", flush=True)
    return row


def write_csv(step, rows):
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, f"step{step}.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("wrote", path, flush=True)


def default_o(L):
    return math.ceil(L / 2)


def snps_for(N, ratio=SNP_RATIO):
    return round(N * ratio)


def step0():
    return [run_case(0, "default", "-", t, N_BASE, L_DEFAULT, COV_DEFAULT,
                     snps_for(N_BASE), default_o(L_DEFAULT)) for t in TOOLS]


def step1():
    rows = []
    for M in M_STEP1:
        N = M * L_DEFAULT // COV_DEFAULT
        rows.append(run_case(1, "M", M, "ssake", N, L_DEFAULT, COV_DEFAULT,
                             snps_for(N), default_o(L_DEFAULT)))
        if M == M_STEP1[0]:
            rows.append(run_case(1, "M", M, "greedy_overlap", N, L_DEFAULT, COV_DEFAULT,
                                 snps_for(N), default_o(L_DEFAULT)))
    return rows


def step2():
    rows = []
    for cov in [5, 10, 15]:
        for t in TOOLS:
            rows.append(run_case(2, "coverage", cov, t, N_BASE, L_DEFAULT, cov,
                                 snps_for(N_BASE), default_o(L_DEFAULT)))
    return rows


def step3():
    rows = []
    for ratio in [0.001, 0.005, 0.01]:
        for t in TOOLS:
            rows.append(run_case(3, "snp_ratio", ratio, t, N_BASE, L_DEFAULT, COV_DEFAULT,
                                 snps_for(N_BASE, ratio), default_o(L_DEFAULT)))
    return rows


def step4():
    rows = []
    for L in [32, 50, 75, 100]:
        for t in TOOLS:
            rows.append(run_case(4, "L", L, t, N_BASE, L, COV_DEFAULT,
                                 snps_for(N_BASE), default_o(L)))
    for o in [20, 25, 35]:
        for t in TOOLS:
            rows.append(run_case(4, "o", o, t, N_BASE, L_DEFAULT, COV_DEFAULT,
                                 snps_for(N_BASE), o))
    return rows


def step5():
    rows = []
    for err in [0, 0.005, 0.01, 0.02]:
        for t in TOOLS:
            rows.append(run_case(5, "err", err, t, N_BASE, L_DEFAULT, COV_DEFAULT,
                                 snps_for(N_BASE), default_o(L_DEFAULT), err=err))
    for d in [0, 1, 2, 3]:
        for t in TOOLS:
            rows.append(run_case(5, "D", d, t, N_BASE, L_DEFAULT, COV_DEFAULT,
                                 snps_for(N_BASE), default_o(L_DEFAULT), D=d, err=0.01))
    return rows


def step6():
    rows = []
    for rl in [25, 50, 100, 200]:
        for t in TOOLS:
            rows.append(run_case(6, "repeat_len", rl, t, N_BASE, L_DEFAULT, COV_DEFAULT,
                                 snps_for(N_BASE), default_o(L_DEFAULT),
                                 repeats=20, repeat_len=rl))
    return rows


STEPS = {0: step0, 1: step1, 2: step2, 3: step3, 4: step4, 5: step5, 6: step6}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", nargs="+", type=int, default=list(range(7)))
    args = parser.parse_args()

    build_tools()
    for s in args.steps:
        write_csv(s, STEPS[s]())


if __name__ == "__main__":
    sys.exit(main())
