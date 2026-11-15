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

// Replaces `count` distinct positions of seq with a different base. Positions
// are drawn by rejection, so the sparse case needs no N-sized array.
inline std::vector<Snp> insert_snps(std::string& seq, int count, std::mt19937_64& rng) {
    const size_t n = seq.size();
    static const char kBase[4] = {'A', 'C', 'G', 'T'};

    std::unordered_set<size_t> chosen;
    std::vector<Snp> snps;
    snps.reserve(count);
    while ((int)snps.size() < count) {
        size_t p = rng() % n;
        if (!chosen.insert(p).second) continue;
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

// Copies ref into a sample genome and places snp_count SNPs in it.
inline Dataset make_sample(const std::string& ref, int L, int snp_count, std::mt19937_64& rng) {
    require_acgt(ref);
    if (L <= 0 || static_cast<size_t>(L) > ref.size())
        throw std::invalid_argument("read length must be in [1, N]");
    if (snp_count < 0 || static_cast<size_t>(snp_count) > ref.size())
        throw std::invalid_argument("snp count must be in [0, N]");
    Dataset d;
    d.sample = ref;
    d.snps = insert_snps(d.sample, snp_count, rng);
    return d;
}

// Samples M = round(N * coverage / L) reads of length L at uniform start positions
// in [0, N - L]. Each base is replaced by a different base with probability err_rate.
// Calls emit(start, read) per read so callers need not hold all reads in memory.
template <class Emit>
void sample_reads(const std::string& sample, int L, double coverage, double err_rate,
                  std::mt19937_64& rng, Emit emit) {
    if (!(err_rate >= 0 && err_rate < 1)) throw std::invalid_argument("err rate must be in [0, 1)");
    if (L <= 0 || static_cast<size_t>(L) > sample.size())
        throw std::invalid_argument("read length must be in [1, N]");
    static const char kBase[4] = {'A', 'C', 'G', 'T'};
    const size_t n = sample.size();
    const size_t m = static_cast<size_t>(std::llround(n * coverage / L));
    std::uniform_int_distribution<size_t> start_dist(0, n - L);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::string read;
    for (size_t i = 0; i < m; ++i) {
        size_t st = start_dist(rng);
        read = sample.substr(st, L);
        if (err_rate > 0) {
            for (char& c : read) {
                if (unit(rng) >= err_rate) continue;
                int b = 0;
                while (kBase[b] != c) ++b;
                c = kBase[(b + 1 + rng() % 3) % 4];
            }
        }
        emit(st, read);
    }
}

// Builds the sample genome and keeps all M reads in memory (tests and small runs).
inline Dataset make_dataset(const std::string& ref, int L, double coverage,
                            int snp_count, uint64_t seed, double err_rate = 0.0) {
    auto rng = make_rng(seed);
    Dataset d = make_sample(ref, L, snp_count, rng);
    sample_reads(d.sample, L, coverage, err_rate, rng,
                 [&](size_t, const std::string& r) { d.reads.push_back(r); });
    return d;
}

#ifndef GEN_READS_NO_MAIN
int main(int argc, char** argv) {
    std::string ref_path, out_dir;
    int L = 30, snps = 50;
    double coverage = 10, err = 0;
    uint64_t seed = 1;
    for (int a = 1; a + 1 < argc; a += 2) {
        std::string flag = argv[a];
        std::string val = argv[a + 1];
        if (flag == "--ref") ref_path = val;
        else if (flag == "--L") L = std::stoi(val);
        else if (flag == "--coverage") coverage = std::stod(val);
        else if (flag == "--snps") snps = std::stoi(val);
        else if (flag == "--err") err = std::stod(val);
        else if (flag == "--seed") seed = std::stoull(val);
        else if (flag == "--out") out_dir = val;
        else { std::cerr << "unknown flag " << flag << "\n"; return 1; }
    }
    if (ref_path.empty() || out_dir.empty()) {
        std::cerr << "usage: gen_reads --ref FILE --L L --coverage C --snps K --err P --seed S --out DIR\n";
        return 1;
    }
    try {
        auto rng = make_rng(seed);
        Dataset d = make_sample(read_fasta_seq(ref_path), L, snps, rng);
        std::filesystem::create_directories(out_dir);

        std::ofstream(out_dir + "/sample.txt") << d.sample << "\n";

        std::ofstream truth(out_dir + "/truth_snps.tsv");
        for (const Snp& s : d.snps) truth << s.pos << '\t' << s.ref << '\t' << s.alt << '\n';

        std::ofstream reads(out_dir + "/reads.txt");
        sample_reads(d.sample, L, coverage, err, rng,
                     [&](size_t, const std::string& r) { reads << r << '\n'; });
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
