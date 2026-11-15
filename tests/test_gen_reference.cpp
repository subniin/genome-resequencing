#define GEN_REFERENCE_NO_MAIN
#include "../src/gen_reference.cpp"
#include <cassert>
#include <unordered_map>
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

    auto rng9 = make_rng(9);
    std::string s = gen_unique_kmer_seq(20000, 16, rng9);
    std::string before = s;
    insert_repeats(s, 0, 100, rng9);
    assert(s == before);
    insert_repeats(s, 3, 100, rng9);
    assert(s.size() == 20000);
    for (char c : s) assert(c == 'A' || c == 'C' || c == 'G' || c == 'T');
    // Each copy makes at least one 100-mer appear twice.
    std::unordered_map<std::string, int> cnt;
    for (size_t i = 0; i + 100 <= s.size(); ++i) cnt[s.substr(i, 100)]++;
    int dup = 0;
    for (auto& kv : cnt) dup += kv.second > 1;
    assert(dup >= 3);
    threw = false;
    try { insert_repeats(s, 1, 15000, rng9); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    return 0;
}
