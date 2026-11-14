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
    return 0;
}
