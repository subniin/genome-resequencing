// Generates a random ACGT reference in which every k-mer is unique.
// Build: g++ -O2 -std=c++17 -static src/gen_reference.cpp -o build/gen_reference.exe
#include "common.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

inline int base_index(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        default:  return 3; // 'T'
    }
}

// Code of the k-mer that ends at index `end` (requires end + 1 >= k).
inline uint64_t kmer_code(const std::string& seq, size_t end, int k) {
    uint64_t code = 0;
    for (size_t t = end + 1 - k; t <= end; ++t) code = (code << 2) | base_index(seq[t]);
    return code;
}

// Builds a length-N sequence where all k-mers are distinct. Bases are placed
// one by one in a random order; when no base fits, the previous base is
// removed (unmarking its k-mer) and tried with its next candidate.
inline std::string gen_unique_kmer_seq(size_t N, int k, std::mt19937_64& rng) {
    if (k < 1 || k > 15) throw std::invalid_argument("k must be in [1, 15]");
    const uint64_t space = 1ULL << (2 * k);
    if (N > space) throw std::invalid_argument("N exceeds 4^k; k-mers cannot be unique");

    static const char kBase[4] = {'A', 'C', 'G', 'T'};
    std::vector<bool> used(space, false);
    std::string seq;
    seq.reserve(N);
    std::vector<uint8_t> off;   // random starting base per position
    std::vector<uint8_t> tried; // how many candidates were tried per position

    while (seq.size() < N) {
        size_t i = seq.size();
        if (off.size() == i) {
            off.push_back((uint8_t)(rng() % 4));
            tried.push_back(0);
        }
        bool placed = false;
        while (tried[i] < 4) {
            int b = (off[i] + tried[i]) % 4;
            tried[i]++;
            seq.push_back(kBase[b]);
            if (i + 1 < (size_t)k) { placed = true; break; }
            uint64_t code = kmer_code(seq, i, k);
            if (!used[code]) {
                used[code] = true;
                placed = true;
                break;
            }
            seq.pop_back();
        }
        if (placed) continue;

        // Dead end at position i: drop it and backtrack to i-1.
        if (i == 0) throw std::runtime_error("no unique-kmer sequence found");
        off.pop_back();
        tried.pop_back();
        size_t j = i - 1;
        if (j + 1 >= (size_t)k) used[kmer_code(seq, j, k)] = false;
        seq.pop_back();
    }
    return seq;
}

#ifndef GEN_REFERENCE_NO_MAIN
int main(int argc, char** argv) {
    size_t N = 10000;
    int k = 14;
    uint64_t seed = 1;
    std::string out_path = "data/ref_10000.txt";
    for (int a = 1; a + 1 < argc; a += 2) {
        std::string flag = argv[a];
        std::string val = argv[a + 1];
        if (flag == "--N") N = std::stoull(val);
        else if (flag == "--k") k = std::stoi(val);
        else if (flag == "--seed") seed = std::stoull(val);
        else if (flag == "--out") out_path = val;
        else { std::cerr << "unknown flag " << flag << "\n"; return 1; }
    }
    try {
        auto rng = make_rng(seed);
        std::string seq = gen_unique_kmer_seq(N, k, rng);
        std::ofstream out(out_path);
        if (!out) { std::cerr << "cannot write " << out_path << "\n"; return 1; }
        out << seq << "\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif
