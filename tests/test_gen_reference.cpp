#define GEN_REFERENCE_NO_MAIN
#include "../src/gen_reference.cpp"
#include <cassert>
#include <unordered_set>

int main() {
    const size_t N = 10000;
    const int k = 14;

    auto rng = make_rng(1);
    std::string seq = gen_unique_kmer_seq(N, k, rng);
    assert(seq.size() == N);
    for (char c : seq) assert(c == 'A' || c == 'C' || c == 'G' || c == 'T');

    std::unordered_set<std::string> kmers;
    for (size_t i = 0; i + k <= seq.size(); ++i) {
        bool inserted = kmers.insert(seq.substr(i, k)).second;
        assert(inserted);
    }

    auto rng2 = make_rng(1);
    assert(gen_unique_kmer_seq(N, k, rng2) == seq);

    auto rng3 = make_rng(1);
    bool threw = false;
    try {
        gen_unique_kmer_seq(5, 1, rng3);
    } catch (const std::exception&) {
        threw = true;
    }
    assert(threw);
    return 0;
}
