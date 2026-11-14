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
    return 0;
}
