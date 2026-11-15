"""Evaluates one fixed set of SSAKE contigs with several evaluator anchor lengths k.
Setup (the old README 4.9 case): N=1e7, L=30, o=15, coverage 10, SNP 0.5%, D=2, no read errors.
Usage: python -u scripts/check_eval_k.py   ->  results/eval_k_check.csv
"""
import csv
import os
import subprocess

import run_experiments as R

N, L, COV, O = 10_000_000, 30, 10, 15
KS = [14, 20, 24, 28]


def main():
    R.build_tools()
    ref = R.ensure_reference(N, 0, 0)
    run_dir = os.path.join(R.DATA, "runs", "eval_k_check")
    os.makedirs(run_dir, exist_ok=True)
    subprocess.run([R.exe("gen_reads"), "--ref", ref, "--L", str(L), "--coverage", str(COV),
                    "--snps", str(R.snps_for(N, 0.005)), "--seed", str(R.SEED), "--out", run_dir],
                   check=True)
    contigs = os.path.join(run_dir, "contigs.fa")
    subprocess.run([R.exe("ssake"), "--reads", os.path.join(run_dir, "reads.txt"), "--o", str(O),
                    "--D", str(R.D), "--out", contigs], check=True)
    rows = []
    for k in KS:
        out = subprocess.run(
            [R.exe("evaluate"), "--contigs", contigs, "--ref", ref,
             "--sample", os.path.join(run_dir, "sample.txt"),
             "--truth", os.path.join(run_dir, "truth_snps.tsv"), "--k", str(k)],
            check=True, capture_output=True, text=True).stdout
        row = {"k": k}
        row.update(dict(line.split("=", 1) for line in out.strip().splitlines()))
        rows.append(row)
        print(row, flush=True)
    path = os.path.join(R.RESULTS, "eval_k_check.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", path)


if __name__ == "__main__":
    main()
