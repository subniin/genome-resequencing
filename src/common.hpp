#pragma once
#include <cctype>
#include <cstdint>
#include <fstream>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <string>
#include <vector>

// Reads the sequence of a FASTA (or plain text) file: skips '>' header lines,
// drops CR/LF, upper-cases bases.
inline std::string read_fasta_seq(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string line, seq;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') continue;
        for (char c : line) {
            if (c == '\r' || c == '\n') continue;
            seq.push_back((char)std::toupper((unsigned char)c));
        }
    }
    return seq;
}

// Reads non-empty lines, stripping trailing CR.
inline std::vector<std::string> read_lines(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<std::string> out;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) out.push_back(line);
    }
    return out;
}

// Writes contigs as FASTA with names contig_0, contig_1, ...
inline void write_fasta(const std::string& path, const std::vector<std::string>& contigs) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    for (size_t i = 0; i < contigs.size(); ++i) {
        out << ">contig_" << i << "\n" << contigs[i] << "\n";
    }
}

// Counts mismatches between a[ai..ai+len) and b[bi..bi+len).
// Stops at limit+1 as soon as the count exceeds limit.
inline int mismatches_bounded(const std::string& a, size_t ai,
                              const std::string& b, size_t bi,
                              size_t len, int limit) {
    int mm = 0;
    for (size_t i = 0; i < len; ++i) {
        if (a[ai + i] != b[bi + i] && ++mm > limit) return limit + 1;
    }
    return mm;
}

inline std::mt19937_64 make_rng(uint64_t seed) {
    return std::mt19937_64(seed);
}

// Keeps the first copy of each distinct read, in input order.
inline std::vector<std::string> dedupe_reads(const std::vector<std::string>& reads) {
    std::unordered_set<std::string> seen;
    std::vector<std::string> out;
    for (const std::string& r : reads) {
        if (seen.insert(r).second) out.push_back(r);
    }
    return out;
}
