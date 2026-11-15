// Evaluates contigs against the reference (SNP calls) and the sample (accuracy).
// Build: g++ -O2 -std=c++17 -static src/evaluate.cpp -o build/evaluate.exe
#include "common.hpp"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

struct TruthSnp {
    size_t pos;
    char alt;
};

struct Metrics {
    double accuracy = 0;
    double snp_recall = 0;
    double snp_precision = 0;
    int n_contigs = 0;
    size_t n50 = 0;
    double covered_frac = 0;
    int misassemblies = 0;
    int unaligned_contigs = 0;
};

// A segment is a run of anchors on one diagonal d = ref_pos - contig_pos.
struct Anchor {
    int cpos;
    int rpos;
};

inline size_t compute_n50(std::vector<size_t> lens) {
    std::sort(lens.rbegin(), lens.rend());
    size_t total = 0;
    for (size_t l : lens) total += l;
    size_t acc = 0;
    for (size_t l : lens) {
        acc += l;
        if (acc * 2 >= total) return l;
    }
    return 0;
}

// Aligns contigs to ref with unique k-mer anchors. Within each diagonal segment,
// contig bases are compared to the reference (SNP calls) and to the sample (accuracy).
inline Metrics evaluate_contigs(const std::vector<std::string>& contigs,
                                const std::string& ref, const std::string& sample,
                                const std::vector<TruthSnp>& truth, int k) {
    if (k < 1 || k > 32) throw std::invalid_argument("k must be in [1, 32]");
    const uint64_t mask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);

    // Unique reference k-mers as (code, position), sorted by code; repeated codes are dropped.
    std::vector<std::pair<uint64_t, uint32_t>> refK;
    if (ref.size() >= (size_t)k) refK.reserve(ref.size() - k + 1);
    uint64_t code = 0;
    for (size_t q = 0; q < ref.size(); ++q) {
        code = ((code << 2) | (uint64_t)base_index(ref[q])) & mask;
        if (q + 1 >= (size_t)k) refK.push_back({code, (uint32_t)(q + 1 - k)});
    }
    std::sort(refK.begin(), refK.end());
    size_t kept = 0;
    for (size_t i = 0; i < refK.size();) {
        size_t j = i + 1;
        while (j < refK.size() && refK[j].first == refK[i].first) ++j;
        if (j - i == 1) refK[kept++] = refK[i];
        i = j;
    }
    refK.resize(kept);

    std::map<size_t, char> truthMap;
    for (const TruthSnp& t : truth) truthMap[t.pos] = t.alt;

    Metrics m;
    m.n_contigs = (int)contigs.size();
    std::vector<size_t> lens;
    size_t totalBases = 0, matchedBases = 0;
    std::vector<bool> covered(ref.size(), false);
    size_t called = 0, hits = 0;
    std::map<size_t, bool> calledSeen; // a SNP is counted once even if two segments cover it

    for (const std::string& c : contigs) {
        lens.push_back(c.size());
        totalBases += c.size();

        std::vector<Anchor> anchors;
        uint64_t cc = 0;
        for (int i = 0; i < (int)c.size(); ++i) {
            cc = ((cc << 2) | (uint64_t)base_index(c[i])) & mask;
            if (i + 1 < k) continue;
            auto it = std::lower_bound(refK.begin(), refK.end(), std::make_pair(cc, (uint32_t)0));
            if (it != refK.end() && it->first == cc) anchors.push_back({i + 1 - k, (int)it->second});
        }
        // Drop isolated spurious anchors: diagonal differs from both neighbours.
        // Ends are kept only if they agree with their single neighbour.
        if (anchors.size() > 1) {
            auto diag = [](const Anchor& x) { return x.rpos - x.cpos; };
            std::vector<Anchor> kept;
            for (size_t i = 0; i < anchors.size(); ++i) {
                bool prevOk = i > 0 && diag(anchors[i]) == diag(anchors[i - 1]);
                bool nextOk = i + 1 < anchors.size() && diag(anchors[i]) == diag(anchors[i + 1]);
                if (prevOk || nextOk) kept.push_back(anchors[i]);
            }
            anchors.swap(kept);
        }
        if (anchors.empty()) {
            m.unaligned_contigs++;
            continue;
        }

        // Split anchors into segments of constant diagonal. A contig base is scored
        // once: a segment never starts before the previous segment's end.
        size_t start = 0;
        int prevEnd = 0;
        for (size_t a = 1; a <= anchors.size(); ++a) {
            bool split = a == anchors.size() ||
                         (anchors[a].rpos - anchors[a].cpos) != (anchors[a - 1].rpos - anchors[a - 1].cpos);
            if (!split) continue;
            if (a < anchors.size()) m.misassemblies++;

            int d = anchors[start].rpos - anchors[start].cpos;
            int cFrom = std::max(start == 0 ? 0 : anchors[start].cpos, prevEnd); // first segment extends to contig start
            int cTo = a == anchors.size() ? (int)c.size() : anchors[a - 1].cpos + k; // exclusive; last segment extends to contig end
            prevEnd = cTo;
            for (int i = cFrom; i < cTo; ++i) {
                int q = i + d;
                if (q < 0 || q >= (int)ref.size()) continue;
                covered[q] = true;
                if (c[i] == sample[q]) matchedBases++;
                if (c[i] != ref[q]) {
                    // Called SNP: reference position q, alternative base c[i].
                    if (!calledSeen.count(q)) {
                        calledSeen[q] = true;
                        called++;
                        auto t = truthMap.find(q);
                        if (t != truthMap.end() && t->second == c[i]) hits++;
                    }
                }
            }
            start = a;
        }
    }

    m.accuracy = totalBases ? (double)matchedBases / totalBases : 0;
    m.snp_recall = truth.empty() ? 0 : (double)hits / truth.size();
    m.snp_precision = called ? (double)hits / called : 0;
    m.n50 = compute_n50(lens);
    size_t cov = std::count(covered.begin(), covered.end(), true);
    m.covered_frac = ref.empty() ? 0 : (double)cov / ref.size();
    return m;
}

#ifndef EVALUATE_NO_MAIN
int main(int argc, char** argv) {
    std::string contigs_path, ref_path, sample_path, truth_path;
    int k = 24;
    for (int a = 1; a + 1 < argc; a += 2) {
        std::string flag = argv[a];
        std::string val = argv[a + 1];
        if (flag == "--contigs") contigs_path = val;
        else if (flag == "--ref") ref_path = val;
        else if (flag == "--sample") sample_path = val;
        else if (flag == "--truth") truth_path = val;
        else if (flag == "--k") k = std::stoi(val);
        else { std::cerr << "unknown flag " << flag << "\n"; return 1; }
    }
    if (contigs_path.empty() || ref_path.empty() || sample_path.empty() || truth_path.empty()) {
        std::cerr << "usage: evaluate --contigs F --ref F --sample F --truth TSV [--k K (default 24)]\n";
        return 1;
    }
    try {
        std::vector<std::string> contigs;
        {
            std::ifstream in(contigs_path);
            if (!in) throw std::runtime_error("cannot open " + contigs_path);
            std::string line, cur;
            bool have = false;
            while (std::getline(in, line)) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (!line.empty() && line[0] == '>') {
                    if (have) contigs.push_back(cur);
                    cur.clear();
                    have = true;
                } else {
                    cur += line;
                }
            }
            if (have) contigs.push_back(cur);
        }
        std::vector<TruthSnp> truth;
        {
            std::ifstream in(truth_path);
            if (!in) throw std::runtime_error("cannot open " + truth_path);
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty()) continue;
                size_t p = std::stoull(line.substr(0, line.find('\t')));
                size_t t2 = line.rfind('\t');
                truth.push_back({p, line[t2 + 1]});
            }
        }
        Metrics m = evaluate_contigs(contigs, read_fasta_seq(ref_path), read_fasta_seq(sample_path), truth, k);
        std::printf("accuracy=%.6f\nsnp_recall=%.6f\nsnp_precision=%.6f\nn_contigs=%d\nn50=%zu\ncovered_frac=%.6f\nmisassemblies=%d\nunaligned_contigs=%d\n",
                    m.accuracy, m.snp_recall, m.snp_precision, m.n_contigs, m.n50,
                    m.covered_frac, m.misassemblies, m.unaligned_contigs);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
