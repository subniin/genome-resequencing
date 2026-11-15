#define SSAKE_NO_MAIN
#include "../src/ssake.cpp"
#include <cassert>
#include <set>

int main() {
    // Same genome and reads as the greedy test: 300 bases, five reads of 100 at shifts of 50.
    auto rng = make_rng(3);
    std::string genome;
    for (int i = 0; i < 300; ++i) genome.push_back("ACGT"[rng() % 4]);
    std::vector<std::string> reads;
    for (int start = 0; start <= 200; start += 50) reads.push_back(genome.substr(start, 100));

    std::vector<std::string> contigs = assemble_ssake(reads, 15, 2);
    assert(contigs.size() == 1);
    assert(contigs[0] == genome);

    // Read order must not change the contig set.
    std::vector<std::string> shuffled = {reads[3], reads[0], reads[4], reads[2], reads[1]};
    std::vector<std::string> c2 = assemble_ssake(shuffled, 15, 2);
    assert(c2.size() == 1);
    assert(c2[0] == genome);

    // Identical copies must not loop forever and must not add contigs.
    std::vector<std::string> dup = {reads[0], reads[1], reads[2], reads[2], reads[2], reads[3], reads[4]};
    std::vector<std::string> c3 = assemble_ssake(dup, 15, 2);
    assert(c3.size() == 1);
    assert(c3[0] == genome);

    // Two reads with no overlap -> two contigs.
    std::vector<std::string> disjoint = {std::string(30, 'A'), std::string(30, 'C')};
    assert(assemble_ssake(disjoint, 15, 2).size() == 2);

    // o >= L is rejected.
    bool threw = false;
    try { assemble_ssake(reads, 100, 2); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    { PackedReads p; p.push("ACGTACGT"); p.push("CGTACGTT"); KeyIndex ix = build_prefix_index(p, 4); auto rg = ix.range(kmer_code("CGTA", 3, 4)); assert(rg.second - rg.first == 1 && ix.ids[rg.first] == 1); }

    // L=100, o=50 (> 32): key is capped, full overlap still checked.
    {
        auto r = make_rng(21);
        std::string g;
        for (int i = 0; i < 1000; ++i) g.push_back("ACGT"[r() % 4]);
        std::vector<std::string> rd;
        for (int st = 0; st + 100 <= 1000; st += 30) rd.push_back(g.substr(st, 100));
        rd.push_back(g.substr(900, 100));
        std::vector<std::string> c = assemble_ssake(rd, 50, 2);
        assert(c.size() == 1 && c[0] == g);
    }
    // Error inside a read's key region: no crash, contigs stay substrings of the genome
    // except at most D mismatches each.
    {
        auto r = make_rng(22);
        std::string g;
        for (int i = 0; i < 400; ++i) g.push_back("ACGT"[r() % 4]);
        std::vector<std::string> rd;
        for (int st = 0; st + 50 <= 400; st += 10) rd.push_back(g.substr(st, 50));
        rd[5][3] = rd[5][3] == 'A' ? 'C' : 'A';
        std::vector<std::string> c = assemble_ssake(rd, 25, 2);
        assert(!c.empty());
    }
    // Empty and single-read inputs.
    assert(assemble_ssake(std::vector<std::string>{}, 15, 2).empty());
    {
        std::vector<std::string> one = {std::string(30, 'A')};
        std::vector<std::string> c = assemble_ssake(one, 15, 2);
        assert(c.size() == 1 && c[0] == one[0]);
    }

    // Coverage 6 reads of a 400-base genome; one read carries an error at its 3' end.
    // vote=true must reproduce the genome exactly.
    {
        auto r = make_rng(41);
        std::string g;
        for (int i = 0; i < 400; ++i) g.push_back("ACGT"[r() % 4]);
        PackedReads p;
        for (int st = 0; st + 50 <= 400; st += 8) {
            std::string rd = g.substr(st, 50);
            if (st == 80) rd[49] = rd[49] == 'A' ? 'C' : 'A';
            p.push(rd);
        }
        p.push(g.substr(350, 50));
        std::vector<std::string> c = assemble_ssake(p, 25, 2, true);
        assert(c.size() == 1 && c[0] == g);
    }
    return 0;
}
