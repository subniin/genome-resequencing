// Benchmark: greedy overlap assembly with all-pairs overlap search.
// Build: g++ -O2 -std=c++17 -static src/greedy_overlap.cpp -o build/greedy_overlap.exe
#include "common.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// Best overlap of a's suffix with b's prefix, length in [o, min(|a|,|b|)-1],
// allowing at most D mismatches. Returns 0 when none.
inline int best_overlap(const std::string& a, const std::string& b, int o, int D) {
    int kmax = (int)std::min(a.size(), b.size()) - 1;
    for (int k = kmax; k >= o; --k) {
        if (mismatches_bounded(a, a.size() - k, b, 0, k, D) <= D) return k;
    }
    return 0;
}

inline int uf_find(std::vector<int>& p, int x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
}

// Greedy assembly: all pairs are scored once, then links are added from the
// longest overlap down. A link is kept when both ends are free and it does
// not close a cycle. Each chain of linked reads is one contig.
inline std::vector<std::string> assemble_greedy(const std::vector<std::string>& input, int o, int D) {
    if (o < 1) throw std::invalid_argument("min overlap must be >= 1");
    for (const std::string& r : input) {
        if (o >= (int)r.size()) throw std::invalid_argument("min overlap must be shorter than the reads");
    }
    const std::vector<std::string> reads = dedupe_reads(input);
    const int m = (int)reads.size();

    // Only pairs with a valid overlap are kept; random data leaves most pairs out.
    std::vector<std::tuple<int, int, int>> cand; // (overlap, i, j)
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i == j) continue;
            int k = best_overlap(reads[i], reads[j], o, D);
            if (k > 0) cand.emplace_back(k, i, j);
        }
    }
    std::sort(cand.begin(), cand.end(), [](const auto& x, const auto& y) {
        if (std::get<0>(x) != std::get<0>(y)) return std::get<0>(x) > std::get<0>(y);
        if (std::get<1>(x) != std::get<1>(y)) return std::get<1>(x) < std::get<1>(y);
        return std::get<2>(x) < std::get<2>(y);
    });

    std::vector<int> next(m, -1), prev(m, -1), ov(m, 0), uf(m);
    for (int i = 0; i < m; ++i) uf[i] = i;
    for (const auto& [k, i, j] : cand) {
        if (next[i] != -1 || prev[j] != -1) continue;
        int ri = uf_find(uf, i), rj = uf_find(uf, j);
        if (ri == rj) continue; // would close a cycle
        next[i] = j;
        prev[j] = i;
        ov[i] = k;
        uf[ri] = rj;
    }

    std::vector<std::string> contigs;
    for (int s = 0; s < m; ++s) {
        if (prev[s] != -1) continue;
        std::string contig = reads[s];
        for (int cur = s; next[cur] != -1; cur = next[cur]) {
            contig += reads[next[cur]].substr(ov[cur]);
        }
        contigs.push_back(contig);
    }
    return contigs;
}

#ifndef GREEDY_OVERLAP_NO_MAIN
int main(int argc, char** argv) {
    std::string reads_path, out_path;
    int o = 15, D = 2;
    for (int a = 1; a + 1 < argc; a += 2) {
        std::string flag = argv[a];
        std::string val = argv[a + 1];
        if (flag == "--reads") reads_path = val;
        else if (flag == "--o") o = std::stoi(val);
        else if (flag == "--D") D = std::stoi(val);
        else if (flag == "--out") out_path = val;
        else { std::cerr << "unknown flag " << flag << "\n"; return 1; }
    }
    if (reads_path.empty() || out_path.empty()) {
        std::cerr << "usage: greedy_overlap --reads FILE --o O --D D --out FASTA\n";
        return 1;
    }
    try {
        write_fasta(out_path, assemble_greedy(read_lines(reads_path), o, D));
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
