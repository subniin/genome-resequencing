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
// vote=false: the first valid read at the largest overlap extends the contig in one step.
// vote=true: every valid unused read (any overlap >= o) votes on the next base, each
// read at its own offset. A base is added with at least 2 voters and a top share of
// at least 2/3, or when a single voter is left (no conflict to resolve).
inline std::vector<std::string> assemble_ssake(PackedReads reads, int o, int D, bool vote) {
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

        // Caps contig growth so a repeat-driven vote cannot loop forever.
        const size_t cap = m * (size_t)L;

        // 3' extension: reads whose prefix overlaps the contig's end.
        for (bool extended = true; extended;) {
            extended = false;
            const int kmax = std::min((int)contig.size(), L - 1);
            std::vector<size_t> cand; // read ids
            std::vector<int> ov;      // overlap length of each candidate with the contig end
            for (int k = kmax; k >= o && (vote || cand.empty()); --k) {
                uint64_t key = kmer_code(contig, contig.size() - k + kk - 1, kk);
                auto [lo, hi] = prefixIdx.range(key);
                for (size_t x = lo; x < hi && (vote || cand.empty()); ++x) {
                    size_t id = prefixIdx.ids[x];
                    if (used[id] || std::find(cand.begin(), cand.end(), id) != cand.end()) continue;
                    if (overlap_ok(reads, id, 0, contig, contig.size() - k, k, D)) {
                        cand.push_back(id);
                        ov.push_back(k);
                    }
                }
            }
            if (cand.empty()) break;

            if (!vote) {
                for (int i = ov[0]; i < L; ++i) contig.push_back(reads.base(cand[0], i));
                used[cand[0]] = true;
                extended = true;
                continue;
            }
            // Vote one base at a time. Candidate i supplies read base ov[i] + t.
            std::vector<bool> active(cand.size(), true);
            int t = 0;
            while (contig.size() < cap) {
                int cnt[4] = {0, 0, 0, 0}, total = 0;
                for (size_t i = 0; i < cand.size(); ++i) {
                    if (active[i] && ov[i] + t < L) {
                        cnt[base_index(reads.base(cand[i], ov[i] + t))]++;
                        total++;
                    }
                }
                if (total == 0) break;
                int best = 0;
                for (int b = 1; b < 4; ++b) if (cnt[b] > cnt[best]) best = b;
                if (total > 1 && cnt[best] * 3 < total * 2) break;
                char bc = "ACGT"[best];
                contig.push_back(bc);
                for (size_t i = 0; i < cand.size(); ++i) {
                    if (active[i] && ov[i] + t < L && reads.base(cand[i], ov[i] + t) != bc) active[i] = false;
                }
                ++t;
            }
            if (t == 0) break;
            extended = true;
            // A candidate is consumed once the vote has passed it: it either finished
            // (its last base is in the contig) or it disagreed on a base and was dropped.
            for (size_t i = 0; i < cand.size(); ++i) {
                if (!active[i] || ov[i] + t >= L) used[cand[i]] = true;
            }
        }

        // 5' extension: reads whose suffix overlaps the contig's start.
        for (bool extended = true; extended;) {
            extended = false;
            const int kmax = std::min((int)contig.size(), L - 1);
            std::vector<size_t> cand;
            std::vector<int> ov;
            for (int k = kmax; k >= o && (vote || cand.empty()); --k) {
                uint64_t key = kmer_code(contig, k - 1, kk);
                auto [lo, hi] = suffixIdx.range(key);
                for (size_t x = lo; x < hi && (vote || cand.empty()); ++x) {
                    size_t id = suffixIdx.ids[x];
                    if (used[id] || std::find(cand.begin(), cand.end(), id) != cand.end()) continue;
                    if (overlap_ok(reads, id, L - k, contig, 0, k, D)) {
                        cand.push_back(id);
                        ov.push_back(k);
                    }
                }
            }
            if (cand.empty()) break;

            if (!vote) {
                std::string pre;
                for (int i = 0; i < L - ov[0]; ++i) pre.push_back(reads.base(cand[0], i));
                contig = pre + contig;
                used[cand[0]] = true;
                extended = true;
                continue;
            }
            // Vote one base at a time, moving left. Candidate i supplies read base
            // L - ov[i] - 1 - t.
            std::vector<bool> active(cand.size(), true);
            int t = 0;
            while (contig.size() < cap) {
                int cnt[4] = {0, 0, 0, 0}, total = 0;
                for (size_t i = 0; i < cand.size(); ++i) {
                    if (active[i] && ov[i] + t < L) {
                        cnt[base_index(reads.base(cand[i], L - ov[i] - 1 - t))]++;
                        total++;
                    }
                }
                if (total == 0) break;
                int best = 0;
                for (int b = 1; b < 4; ++b) if (cnt[b] > cnt[best]) best = b;
                if (total > 1 && cnt[best] * 3 < total * 2) break;
                char bc = "ACGT"[best];
                contig.insert(contig.begin(), bc);
                for (size_t i = 0; i < cand.size(); ++i) {
                    if (active[i] && ov[i] + t < L && reads.base(cand[i], L - ov[i] - 1 - t) != bc) active[i] = false;
                }
                ++t;
            }
            if (t == 0) break;
            extended = true;
            for (size_t i = 0; i < cand.size(); ++i) {
                if (!active[i] || ov[i] + t >= L) used[cand[i]] = true;
            }
        }
        contigs.push_back(std::move(contig));
    }
    return contigs;
}

inline std::vector<std::string> assemble_ssake(PackedReads reads, int o, int D) {
    return assemble_ssake(std::move(reads), o, D, false);
}

// Convenience overload for in-memory reads of one length.
inline std::vector<std::string> assemble_ssake(const std::vector<std::string>& input, int o, int D) {
    PackedReads p;
    for (const std::string& r : input) p.push(r);
    return assemble_ssake(std::move(p), o, D, false);
}

#ifndef SSAKE_NO_MAIN
int main(int argc, char** argv) {
    std::string reads_path, out_path;
    int o = 15, D = 2, vote = 0;
    for (int a = 1; a + 1 < argc; a += 2) {
        std::string flag = argv[a];
        std::string val = argv[a + 1];
        if (flag == "--reads") reads_path = val;
        else if (flag == "--o") o = std::stoi(val);
        else if (flag == "--D") D = std::stoi(val);
        else if (flag == "--vote") vote = std::stoi(val);
        else if (flag == "--out") out_path = val;
        else { std::cerr << "unknown flag " << flag << "\n"; return 1; }
    }
    if (reads_path.empty() || out_path.empty()) {
        std::cerr << "usage: ssake --reads FILE --o O --D D [--vote 0|1] --out FASTA\n";
        return 1;
    }
    try {
        write_fasta(out_path, assemble_ssake(load_reads_packed(reads_path), o, D, vote != 0));
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
