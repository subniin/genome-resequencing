// SSAKE-style seed-and-extend assembler: prefix/suffix hash index, extension in both directions.
// Build: g++ -O2 -std=c++17 -static src/ssake.cpp -o build/ssake.exe
#include "common.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Assembles reads with seed-and-extend. Seeds are taken in input order; each
// seed is extended at its 3' end, then at its 5' end, using only unused reads.
inline std::vector<std::string> assemble_ssake(const std::vector<std::string>& input, int o, int D) {
    if (o < 1) throw std::invalid_argument("min overlap must be >= 1");
    for (const std::string& r : input) {
        if (o >= (int)r.size()) throw std::invalid_argument("min overlap must be shorter than the reads");
    }
    const std::vector<std::string> reads = dedupe_reads(input);
    const int m = (int)reads.size();
    int maxLen = 0;
    for (const std::string& r : reads) maxLen = std::max(maxLen, (int)r.size());

    // prefix index: first o bases -> read ids; suffix index: last o bases -> read ids.
    std::unordered_map<std::string, std::vector<int>> prefixIdx, suffixIdx;
    for (int i = 0; i < m; ++i) {
        prefixIdx[reads[i].substr(0, o)].push_back(i);
        suffixIdx[reads[i].substr(reads[i].size() - o, o)].push_back(i);
    }

    std::vector<bool> used(m, false);
    std::vector<std::string> contigs;

    for (int s = 0; s < m; ++s) {
        if (used[s]) continue;
        used[s] = true;
        std::string contig = reads[s];

        // 3' extension: a read whose prefix overlaps the contig's end.
        for (bool extended = true; extended;) {
            extended = false;
            int kmax = std::min((int)contig.size(), maxLen - 1);
            for (int k = kmax; k >= o && !extended; --k) {
                auto it = prefixIdx.find(contig.substr(contig.size() - k, o));
                if (it == prefixIdx.end()) continue;
                for (int id : it->second) {
                    if (used[id] || (int)reads[id].size() <= k) continue;
                    if (mismatches_bounded(contig, contig.size() - k, reads[id], 0, k, D) <= D) {
                        contig += reads[id].substr(k);
                        used[id] = true;
                        extended = true;
                        break;
                    }
                }
            }
        }

        // 5' extension: a read whose suffix overlaps the contig's start.
        for (bool extended = true; extended;) {
            extended = false;
            int kmax = std::min((int)contig.size(), maxLen - 1);
            for (int k = kmax; k >= o && !extended; --k) {
                auto it = suffixIdx.find(contig.substr(k - o, o));
                if (it == suffixIdx.end()) continue;
                for (int id : it->second) {
                    if (used[id] || (int)reads[id].size() <= k) continue;
                    size_t rs = reads[id].size() - k;
                    if (mismatches_bounded(reads[id], rs, contig, 0, k, D) <= D) {
                        contig = reads[id].substr(0, rs) + contig;
                        used[id] = true;
                        extended = true;
                        break;
                    }
                }
            }
        }
        contigs.push_back(contig);
    }
    return contigs;
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
        write_fasta(out_path, assemble_ssake(read_lines(reads_path), o, D));
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
