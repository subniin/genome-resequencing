#define EVALUATE_NO_MAIN
#include "../src/evaluate.cpp"
#include <cassert>
#include <cmath>

int main() {
    // Reference of 1000 random bases; sample = reference with SNPs at 100, 500, 900.
    auto rng = make_rng(11);
    std::string ref;
    for (int i = 0; i < 1000; ++i) ref.push_back("ACGT"[rng() % 4]);
    std::string sample = ref;
    std::vector<TruthSnp> truth;
    for (size_t p : {100, 500, 900}) {
        char alt = sample[p] == 'A' ? 'C' : 'A';
        sample[p] = alt;
        truth.push_back({p, alt});
    }

    // Contig equal to the sample: perfect recall, precision and accuracy.
    Metrics m = evaluate_contigs({sample}, ref, sample, truth, 14);
    assert(std::fabs(m.snp_recall - 1.0) < 1e-9);
    assert(std::fabs(m.snp_precision - 1.0) < 1e-9);
    assert(std::fabs(m.accuracy - 1.0) < 1e-9);
    assert(m.n_contigs == 1);
    assert(m.misassemblies == 0);
    assert(m.unaligned_contigs == 0);

    // Two halves joined in the wrong order: a misassembly, still aligned.
    std::string swapped = sample.substr(500) + sample.substr(0, 500);
    Metrics w = evaluate_contigs({swapped}, ref, sample, truth, 14);
    assert(w.misassemblies >= 1);
    assert(w.unaligned_contigs == 0);

    // Random chimeras on a large reference: a junction k-mer can match elsewhere,
    // splitting segments. Each contig base must count once, so accuracy stays <= 1.
    auto rng2 = make_rng(5);
    std::string big;
    for (int i = 0; i < 1000000; ++i) big.push_back("ACGT"[rng2() % 4]);
    std::vector<std::string> chimeras;
    for (int t = 0; t < 3000; ++t) {
        size_t a = rng2() % (big.size() - 300);
        size_t b = rng2() % (big.size() - 300);
        chimeras.push_back(big.substr(a, 150) + big.substr(b, 150));
    }
    Metrics c = evaluate_contigs(chimeras, big, big, {}, 14);
    assert(c.accuracy <= 1.0 + 1e-12);

    // Unrelated contig: unaligned, no crash.
    std::string junk(200, 'A');
    Metrics u = evaluate_contigs({junk}, ref, sample, truth, 14);
    assert(u.unaligned_contigs == 1);
    assert(u.snp_recall == 0.0);
    return 0;
}
