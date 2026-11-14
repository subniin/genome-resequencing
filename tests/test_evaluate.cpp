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

    // F1: a SNP-containing 14-mer that also occurs elsewhere in the reference gives an
    // isolated spurious anchor; it must not split the contig or move scoring off-diagonal.
    {
        auto r3 = make_rng(77);
        std::string ref3;
        for (int i = 0; i < 1000000; ++i) ref3.push_back("ACGT"[r3() % 4]);
        std::string smp3 = ref3;
        std::vector<TruthSnp> tr3;
        for (int s = 0; s < 50; ++s) {
            size_t p = 10000 + 19000 * (size_t)s;
            char alt = smp3[p] == 'A' ? 'C' : 'A';
            smp3[p] = alt;
            tr3.push_back({p, alt});
        }
        // Plant the SNP 14-mer of every 5th SNP at a far-away spot (in ref and sample).
        for (int s = 0; s < 50; s += 5) {
            size_t p = tr3[s].pos, x = 500000 + 3000 * (size_t)s;
            for (int j = 0; j < 14; ++j) ref3[x + j] = smp3[x + j] = smp3[p - 6 + j];
            // Flanks differ from the SNP's flanks so the spurious match is a single k-mer.
            ref3[x - 1] = smp3[x - 1] = smp3[p - 7] == 'A' ? 'C' : 'A';
            ref3[x + 14] = smp3[x + 14] = smp3[p + 8] == 'A' ? 'C' : 'A';
        }
        Metrics f1 = evaluate_contigs({smp3}, ref3, smp3, tr3, 14);
        assert(f1.misassemblies == 0);
        assert(std::fabs(f1.snp_recall - 1.0) < 1e-9);
    }

    // F2: SNPs inside the first and last k bases must still be scored (tail bases).
    {
        auto r4 = make_rng(3);
        std::string ref4;
        for (int i = 0; i < 10000; ++i) ref4.push_back("ACGT"[r4() % 4]);
        std::string smp4 = ref4;
        std::vector<TruthSnp> tr4;
        for (size_t p : {3, 5000, 9996}) {
            char alt = smp4[p] == 'A' ? 'C' : 'A';
            smp4[p] = alt;
            tr4.push_back({p, alt});
        }
        Metrics f2 = evaluate_contigs({smp4}, ref4, smp4, tr4, 14);
        assert(f2.accuracy == 1.0);
        assert(std::fabs(f2.snp_recall - 1.0) < 1e-9);
    }

    // Unrelated contig: unaligned, no crash.
    std::string junk(200, 'A');
    Metrics u = evaluate_contigs({junk}, ref, sample, truth, 14);
    assert(u.unaligned_contigs == 1);
    assert(u.snp_recall == 0.0);
    return 0;
}
