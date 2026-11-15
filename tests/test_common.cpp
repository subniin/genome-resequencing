#include "../src/common.hpp"
#include <cassert>
#include <fstream>
#include <cstdio>

int main() {
    assert(mismatches_bounded("ACGT", 0, "ACCT", 0, 4, 2) == 1);
    assert(mismatches_bounded("AAAA", 0, "TTTT", 0, 4, 2) == 3);

    const char* path = "build/tmp_fasta_test.fa";
    {
        std::ofstream out(path, std::ios::binary);
        out << ">hdr\r\nacgt\r\nGG\r\n";
    }
    assert(read_fasta_seq(path) == "ACGTGG");
    std::remove(path);

    // PackedReads round trip and k-mer codes match the string versions.
    PackedReads p;
    std::string a = "ACGTACGTTTGCAACGTACGTTTGCAACGTACGTTTGCAACGTACGTTTGCAACGTACGTTTGCAACA"; // 68 bases, crosses a word
    std::string b = std::string(68, 'G');
    p.push(a); p.push(b); p.push(a);
    assert(p.size() == 3 && p.L == 68);
    assert(p.str(0) == a && p.str(1) == b);
    assert(p.base(0, 67) == a[67]);
    assert(p.code(0, 30, 32) == kmer_code(a, 61, 32));
    bool threw = false;
    try { p.push("ACGT"); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    dedupe_packed(p);
    assert(p.size() == 2 && p.str(0) == a && p.str(1) == b);

    // load_reads_packed: skips blank lines, strips CR, rejects mixed lengths.
    { std::ofstream f("build/tmp_reads.txt", std::ios::binary); f << "acgt\r\n\r\nTTGG\n"; }
    PackedReads q = load_reads_packed("build/tmp_reads.txt");
    assert(q.size() == 2 && q.str(0) == "ACGT" && q.str(1) == "TTGG");
    { std::ofstream f("build/tmp_reads.txt", std::ios::binary); f << "ACGT\nACG\n"; }
    threw = false;
    try { load_reads_packed("build/tmp_reads.txt"); } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);
    std::remove("build/tmp_reads.txt");
    return 0;
}
