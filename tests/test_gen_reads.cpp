#define GEN_READS_NO_MAIN
#include "../src/gen_reads.cpp"
#include <cassert>
#include <cstdio>
#include <fstream>
#include <set>

int main() {
    // insert_snps: 50 distinct positions, each alt differs from ref, sample matches.
    std::string ref;
    for (int i = 0; i < 10000; ++i) ref.push_back("ACGT"[i % 4]);
    std::string sample = ref;
    auto rng = make_rng(7);
    std::vector<Snp> snps = insert_snps(sample, 50, rng);
    assert(snps.size() == 50);
    std::set<size_t> positions;
    for (const Snp& s : snps) {
        assert(s.ref != s.alt);
        assert(ref[s.pos] == s.ref);
        assert(sample[s.pos] == s.alt);
        positions.insert(s.pos);
    }
    assert(positions.size() == 50);
    int diff = 0;
    for (size_t i = 0; i < ref.size(); ++i) diff += ref[i] != sample[i];
    assert(diff == 50);

    // read_fasta_seq normalizes lowercase and CRLF; require_acgt rejects other bases.
    const char* path = "build/tmp_ref_test.fa";
    {
        std::ofstream out(path, std::ios::binary);
        out << ">ref\r\nacgt\r\nAC\r\n";
    }
    std::string norm = read_fasta_seq(path);
    assert(norm == "ACGTAC");
    require_acgt(norm);
    bool threw = false;
    try { require_acgt("ACGNT"); } catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    std::remove(path);

    // make_dataset: deterministic for a fixed seed; default config gives 3333 reads of length 30.
    Dataset a = make_dataset(ref, 30, 10, 50, 1);
    Dataset b = make_dataset(ref, 30, 10, 50, 1);
    assert(a.sample == b.sample);
    assert(a.reads == b.reads);
    assert(a.reads.size() == 3333);
    for (const std::string& r : a.reads) assert(r.size() == 30);
    // F3: SNP count must fit in the reference (and be non-negative).
    {
        std::string tiny(14, 'A');
        bool thr = false;
        try { make_dataset(tiny, 14, 1, 15, 1); } catch (const std::invalid_argument&) { thr = true; }
        assert(thr);
        thr = false;
        try { make_dataset(tiny, 14, 1, -1, 1); } catch (const std::invalid_argument&) { thr = true; }
        assert(thr);
        make_dataset(tiny, 14, 1, 5, 1);
    }

    // err_rate 0: every read equals the sample at its start.
    {
        auto r = make_rng(4);
        std::string sample = ref; // 10000 bases from the existing test
        size_t n = 0;
        sample_reads(sample, 50, 10, 0.0, r, [&](size_t st, const std::string& rd) {
            assert(rd == sample.substr(st, 50)); ++n; });
        assert(n == 2000);
    }
    // err_rate 0.05: mismatch fraction within [0.04, 0.06].
    {
        auto r = make_rng(4);
        size_t mm = 0, total = 0;
        sample_reads(ref, 50, 10, 0.05, r, [&](size_t st, const std::string& rd) {
            for (int i = 0; i < 50; ++i) mm += rd[i] != ref[st + i];
            total += 50; });
        double f = (double)mm / total;
        assert(f > 0.04 && f < 0.06);
    }
    bool thr = false;
    try { make_dataset(ref, 50, 10, 0, 1, 1.5); } catch (const std::invalid_argument&) { thr = true; }
    assert(thr);
    return 0;
}
