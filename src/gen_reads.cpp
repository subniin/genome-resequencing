// Generates a sample genome from a reference (SNPs only) and samples short reads from it.
// Build: g++ -O2 -std=c++17 -static src/gen_reads.cpp -o build/gen_reads.exe
#include "common.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct Snp {
    size_t pos;
    char ref;
    char alt;
};

struct Dataset {
    std::string sample;
    std::vector<Snp> snps;
    std::vector<std::string> reads;
};

// Throws if seq contains anything other than A, C, G, T.
inline void require_acgt(const std::string& seq) {
    for (char c : seq) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            throw std::runtime_error(std::string("non-ACGT base in reference: ") + c);
    }
}

// Replaces `count` distinct positions of seq with a different base.
// Positions are drawn without replacement, so they never collide.
inline std::vector<Snp> insert_snps(std::string& seq, int count, std::mt19937_64& rng) {
    const size_t n = seq.size();
    std::vector<size_t> pos(n);
    for (size_t i = 0; i < n; ++i) pos[i] = i;
    static const char kBase[4] = {'A', 'C', 'G', 'T'};

    std::vector<Snp> snps;
    snps.reserve(count);
    for (int i = 0; i < count; ++i) {
        size_t j = i + rng() % (n - i);
        std::swap(pos[i], pos[j]);
        size_t p = pos[i];
        char ref = seq[p];
        char alt;
        do {
            alt = kBase[rng() % 4];
        } while (alt == ref);
        seq[p] = alt;
        snps.push_back({p, ref, alt});
    }
    return snps;
}

// Builds the sample genome from ref, then samples M = round(N * coverage / L)
// reads of length L at uniform start positions in [0, N - L].
inline Dataset make_dataset(const std::string& ref, int L, double coverage,
                            int snp_count, uint64_t seed) {
    require_acgt(ref);
    if (L <= 0 || static_cast<size_t>(L) > ref.size())
        throw std::invalid_argument("read length must be in [1, N]");
    if (snp_count < 0 || static_cast<size_t>(snp_count) > ref.size())
        throw std::invalid_argument("snp count must be in [0, N]");

    auto rng = make_rng(seed);
    Dataset d;
    d.sample = ref;
    d.snps = insert_snps(d.sample, snp_count, rng);

    const size_t n = d.sample.size();
    const size_t m = static_cast<size_t>(std::llround(n * coverage / L));
    std::uniform_int_distribution<size_t> start_dist(0, n - L);
    d.reads.reserve(m);
    for (size_t i = 0; i < m; ++i) {
        d.reads.push_back(d.sample.substr(start_dist(rng), L));
    }
    return d;
}

#ifndef GEN_READS_NO_MAIN
int main(int argc, char** argv) {
    std::string ref_path, out_dir;
    int L = 30, snps = 50;
    double coverage = 10;
    uint64_t seed = 1;
    for (int a = 1; a + 1 < argc; a += 2) {
        std::string flag = argv[a];
        std::string val = argv[a + 1];
        if (flag == "--ref") ref_path = val;
        else if (flag == "--L") L = std::stoi(val);
        else if (flag == "--coverage") coverage = std::stod(val);
        else if (flag == "--snps") snps = std::stoi(val);
        else if (flag == "--seed") seed = std::stoull(val);
        else if (flag == "--out") out_dir = val;
        else { std::cerr << "unknown flag " << flag << "\n"; return 1; }
    }
    if (ref_path.empty() || out_dir.empty()) {
        std::cerr << "usage: gen_reads --ref FILE --L L --coverage C --snps K --seed S --out DIR\n";
        return 1;
    }
    try {
        Dataset d = make_dataset(read_fasta_seq(ref_path), L, coverage, snps, seed);
        std::filesystem::create_directories(out_dir);

        std::ofstream(out_dir + "/sample.txt") << d.sample << "\n";

        std::ofstream truth(out_dir + "/truth_snps.tsv");
        for (const Snp& s : d.snps) truth << s.pos << '\t' << s.ref << '\t' << s.alt << '\n';

        std::ofstream reads(out_dir + "/reads.txt");
        for (const std::string& r : d.reads) reads << r << '\n';
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
