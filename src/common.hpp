#pragma once
#include <algorithm>
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

inline int base_index(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        default:  return 3; // 'T'
    }
}

// Code of the k-mer that ends at index `end` (requires end + 1 >= k, k in [1, 32]).
inline uint64_t kmer_code(const std::string& seq, size_t end, int k) {
    uint64_t code = 0;
    for (size_t t = end + 1 - k; t <= end; ++t) code = (code << 2) | base_index(seq[t]);
    return code;
}

// Reads stored at 2 bits per base, all reads of one length L. Read i occupies
// `words()` consecutive uint64 words; base j sits at bit 2j of its read.
struct PackedReads {
    int L = 0;
    std::vector<uint64_t> data;
    size_t count = 0;

    size_t size() const { return count; }
    size_t words() const { return (2 * (size_t)L + 63) / 64; }

    // The first read fixes L. Other lengths throw std::invalid_argument;
    // non-ACGT bases throw std::runtime_error.
    void push(const std::string& r) {
        if (r.empty()) throw std::invalid_argument("empty read");
        if (L == 0) L = (int)r.size();
        if ((int)r.size() != L) throw std::invalid_argument("read length differs from L");
        for (char c : r) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                throw std::runtime_error(std::string("non-ACGT base: ") + c);
        }
        const size_t W = words();
        const size_t off = count * W;
        data.resize(off + W, 0);
        for (int i = 0; i < L; ++i) {
            size_t bit = 2 * (size_t)i;
            data[off + bit / 64] |= (uint64_t)base_index(r[i]) << (bit % 64);
        }
        ++count;
    }

    char base(size_t id, int i) const {
        size_t bit = 2 * (size_t)i;
        uint64_t b = (data[id * words() + bit / 64] >> (bit % 64)) & 3;
        return "ACGT"[b];
    }

    // k-mer starting at base `start` of read id, first base in the high bits (same as kmer_code).
    uint64_t code(size_t id, int start, int k) const {
        uint64_t c = 0;
        for (int i = start; i < start + k; ++i) c = (c << 2) | (uint64_t)base_index(base(id, i));
        return c;
    }

    std::string str(size_t id) const {
        std::string s(L, 'A');
        for (int i = 0; i < L; ++i) s[i] = base(id, i);
        return s;
    }
};

// Streams a read file (one read per line) into PackedReads. Empty lines are
// skipped, a trailing CR is removed, bases are upper-cased. A read whose length
// differs from the first one stops loading with the line number in the message.
inline PackedReads load_reads_packed(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    PackedReads p;
    std::string line;
    size_t lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        for (char& c : line) c = (char)std::toupper((unsigned char)c);
        if (p.count && line.size() != (size_t)p.L)
            throw std::invalid_argument("line " + std::to_string(lineNo) + ": read length " +
                                        std::to_string(line.size()) + " differs from " +
                                        std::to_string(p.L));
        try {
            p.push(line);
        } catch (const std::runtime_error& e) {
            throw std::runtime_error("line " + std::to_string(lineNo) + ": " + e.what());
        }
    }
    return p;
}

// Keeps the first copy of each distinct read, in input order. Reads are sorted
// by content (ties by id) to find duplicate groups; only the smallest id of each
// group survives. Temporary memory is one id per read.
inline void dedupe_packed(PackedReads& r) {
    const size_t n = r.count;
    const size_t W = r.words();
    if (n == 0) return;
    std::vector<uint32_t> order(n);
    for (size_t i = 0; i < n; ++i) order[i] = (uint32_t)i;
    auto words_less = [&](uint32_t a, uint32_t b) {
        for (size_t w = 0; w < W; ++w) {
            uint64_t x = r.data[(size_t)a * W + w], y = r.data[(size_t)b * W + w];
            if (x != y) return x < y;
        }
        return a < b;
    };
    std::sort(order.begin(), order.end(), words_less);

    std::vector<bool> keep(n, false);
    for (size_t i = 0; i < n; ++i) {
        bool dup_of_prev = i > 0 && std::equal(r.data.begin() + (size_t)order[i] * W,
                                               r.data.begin() + ((size_t)order[i] + 1) * W,
                                               r.data.begin() + (size_t)order[i - 1] * W);
        if (!dup_of_prev) keep[order[i]] = true;
    }

    size_t out = 0;
    for (size_t i = 0; i < n; ++i) {
        if (!keep[i]) continue;
        if (out != i) std::copy(r.data.begin() + i * W, r.data.begin() + (i + 1) * W, r.data.begin() + out * W);
        ++out;
    }
    r.data.resize(out * W);
    r.count = out;
}
