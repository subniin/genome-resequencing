#define GREEDY_OVERLAP_NO_MAIN
#include "../src/greedy_overlap.cpp"
#include <cassert>

int main() {
    // Genome of 300 random bases; five reads of length 100 at shifts of 50.
    auto rng = make_rng(3);
    std::string genome;
    for (int i = 0; i < 300; ++i) genome.push_back("ACGT"[rng() % 4]);
    std::vector<std::string> reads;
    for (int start = 0; start <= 200; start += 50) reads.push_back(genome.substr(start, 100));

    std::vector<std::string> contigs = assemble_greedy(reads, 15, 2);
    assert(contigs.size() == 1);
    assert(contigs[0] == genome);

    // Read order must not change the contig set.
    std::vector<std::string> shuffled = {reads[3], reads[0], reads[4], reads[2], reads[1]};
    std::vector<std::string> c2 = assemble_greedy(shuffled, 15, 2);
    assert(c2.size() == 1);
    assert(c2[0] == genome);

    // Single read -> single contig equal to that read.
    std::vector<std::string> one = {reads[0]};
    std::vector<std::string> c1 = assemble_greedy(one, 15, 2);
    assert(c1.size() == 1 && c1[0] == reads[0]);

    // Two reads with no overlap -> two contigs.
    std::vector<std::string> disjoint = {std::string(30, 'A'), std::string(30, 'C')};
    assert(assemble_greedy(disjoint, 15, 2).size() == 2);

    // Identical reads terminate and are not merged into a loop.
    std::vector<std::string> dup = {reads[0], reads[0], reads[0]};
    assert(assemble_greedy(dup, 15, 2).size() == 3);

    // o >= L is rejected.
    bool threw = false;
    try { assemble_greedy(reads, 100, 2); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    return 0;
}
