// SSAKE-style seed-and-extend assembler: sorted k-mer key indexes over 2-bit packed reads,
// extension in both directions.
// Build: g++ -O2 -std=c++17 -static src/ssake.cpp -o build/ssake.exe
#include "common.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Keys in ascending order with the read id that carries each one.
// range(key) is the block of entries whose key equals `key`.
struct KeyIndex {
    std::vector<uint64_t> keys;
    std::vector<uint32_t> ids;

    std::pair<size_t, size_t> range(uint64_t key) const {
        auto r = std::equal_range(keys.begin(), keys.end(), key);
        return {size_t(r.first - keys.begin()), size_t(r.second - keys.begin())};
    }
};

// Indexes the first kk bases (suffix == false) or the last kk bases (suffix == true) of each read.
inline KeyIndex build_key_index(const PackedReads& r, int kk, bool suffix) {
    const size_t n = r.size();
    const int start = suffix ? r.L - kk : 0;
    std::vector<std::pair<uint64_t, uint32_t>> pairs(n);
    for (size_t i = 0; i < n; ++i) pairs[i] = {r.code(i, start, kk), (uint32_t)i};
    std::sort(pairs.begin(), pairs.end());
    KeyIndex ix;
    ix.keys.resize(n);
    ix.ids.resize(n);
    for (size_t i = 0; i < n; ++i) {
        ix.keys[i] = pairs[i].first;
        ix.ids[i] = pairs[i].second;
    }
    return ix;
}

inline KeyIndex build_prefix_index(const PackedReads& r, int kk) { return build_key_index(r, kk, false); }
inline KeyIndex build_suffix_index(const PackedReads& r, int kk) { return build_key_index(r, kk, true); }

// True when read `id` bases [rs, rs+len) differ from s[ss, ss+len) in at most D positions.
inline bool overlap_ok(const PackedReads& r, size_t id, int rs, const std::string& s, size_t ss, int len, int D) {
    int mm = 0;
    for (int i = 0; i < len; ++i) {
        if (r.base(id, rs + i) != s[ss + i] && ++mm > D) return false;
    }
    return true;
}

// Assembles reads with seed-and-extend. Seeds are taken in input order; each
// seed is extended at its 3' end, then at its 5' end, using only unused reads.
// The index key is capped at 32 bases; the whole overlap is checked for mismatches.
inline std::vector<std::string> assemble_ssake(PackedReads reads, int o, int D) {
    if (o < 1) throw std::invalid_argument("min overlap must be >= 1");
    if (reads.size() == 0) return {};
    if (o >= reads.L) throw std::invalid_argument("min overlap must be shorter than the reads");
    dedupe_packed(reads);
    const size_t m = reads.size();
    if (m >= UINT32_MAX) throw std::length_error("too many reads for 32-bit read ids");
    const int L = reads.L;
    const int kk = std::min(o, 32);

    KeyIndex prefixIdx = build_prefix_index(reads, kk);
    KeyIndex suffixIdx = build_suffix_index(reads, kk);

    std::vector<bool> used(m, false);
    std::vector<std::string> contigs;

    for (size_t s = 0; s < m; ++s) {
        if (used[s]) continue;
        used[s] = true;
        std::string contig = reads.str(s);

        // 3' extension: a read whose prefix overlaps the contig's end.
        for (bool extended = true; extended;) {
            extended = false;
            int kmax = std::min((int)contig.size(), L - 1);
            for (int k = kmax; k >= o && !extended; --k) {
                uint64_t key = kmer_code(contig, contig.size() - k + kk - 1, kk);
                auto [lo, hi] = prefixIdx.range(key);
                for (size_t x = lo; x < hi && !extended; ++x) {
                    size_t id = prefixIdx.ids[x];
                    if (used[id]) continue;
                    if (overlap_ok(reads, id, 0, contig, contig.size() - k, k, D)) {
                        for (int i = k; i < L; ++i) contig.push_back(reads.base(id, i));
                        used[id] = true;
                        extended = true;
                    }
                }
            }
        }

        // 5' extension: a read whose suffix overlaps the contig's start.
        for (bool extended = true; extended;) {
            extended = false;
            int kmax = std::min((int)contig.size(), L - 1);
            for (int k = kmax; k >= o && !extended; --k) {
                uint64_t key = kmer_code(contig, k - 1, kk);
                auto [lo, hi] = suffixIdx.range(key);
                for (size_t x = lo; x < hi && !extended; ++x) {
                    size_t id = suffixIdx.ids[x];
                    if (used[id]) continue;
                    if (overlap_ok(reads, id, L - k, contig, 0, k, D)) {
                        std::string pre;
                        for (int i = 0; i < L - k; ++i) pre.push_back(reads.base(id, i));
                        contig = pre + contig;
                        used[id] = true;
                        extended = true;
                    }
                }
            }
        }
        contigs.push_back(std::move(contig));
    }
    return contigs;
}

// Convenience overload for in-memory reads of one length.
inline std::vector<std::string> assemble_ssake(const std::vector<std::string>& input, int o, int D) {
    PackedReads p;
    for (const std::string& r : input) p.push(r);
    return assemble_ssake(std::move(p), o, D);
}

#ifndef SSAKE_NO_MAIN
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
        std::cerr << "usage: ssake --reads FILE --o O --D D --out FASTA\n";
        return 1;
    }
    try {
        write_fasta(out_path, assemble_ssake(load_reads_packed(reads_path), o, D));
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
