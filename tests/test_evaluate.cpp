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

    // Unrelated contig: unaligned, no crash.
    std::string junk(200, 'A');
    Metrics u = evaluate_contigs({junk}, ref, sample, truth, 14);
    assert(u.unaligned_contigs == 1);
    assert(u.snp_recall == 0.0);
    return 0;
}
