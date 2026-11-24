// Generates a random ACGT reference in which every k-mer is unique.
// Build: gcc -O2 -std=c11 c/gen_reference.c -o build/c/gen_reference.exe
#include "common.h"

// Builds a length-N sequence where all k-mers are distinct. Bases are placed one by one in
// a random order; when no base fits, the previous base is removed (unmarking its k-mer)
// and tried with its next candidate.
static Str gen_unique_kmer_seq(size_t N, int k, Rng *rng) {
    if (k < 1 || k > 16) die("k must be in [1, 16]");
    const uint64_t space = 1ULL << (2 * k);
    if (N > space) die("N exceeds 4^k; k-mers cannot be unique");

    static const char kBase[4] = {'A', 'C', 'G', 'T'};
    uint64_t *used = calloc(space / 64 + 1, sizeof(uint64_t)); // bitset over k-mer codes
    uint8_t *off = xrealloc(NULL, N + 1), *tried = xrealloc(NULL, N + 1);
    Str seq = {0};
    str_push(&seq, 'A'); seq.n = 0; // allocate

    size_t have = 0; // positions with off/tried set
    while (seq.n < N) {
        size_t i = seq.n;
        if (have == i) {
            off[i] = (uint8_t)(rng_next(rng) % 4);
            tried[i] = 0;
            have++;
        }
        int placed = 0;
        while (tried[i] < 4) {
            int b = (off[i] + tried[i]) % 4;
            tried[i]++;
            str_push(&seq, kBase[b]);
            if (i + 1 < (size_t)k) { placed = 1; break; }
            uint64_t code = kmer_code(seq.s, i, k);
            if (!(used[code / 64] >> (code % 64) & 1)) {
                used[code / 64] |= 1ULL << (code % 64);
                placed = 1;
                break;
            }
            seq.s[--seq.n] = 0;
        }
        if (placed) continue;

        // Dead end at position i: drop it and backtrack to i-1.
        if (i == 0) die("no unique-kmer sequence found");
        have--;
        size_t j = i - 1;
        if (j + 1 >= (size_t)k) {
            uint64_t code = kmer_code(seq.s, j, k);
            used[code / 64] &= ~(1ULL << (code % 64));
        }
        seq.s[--seq.n] = 0;
    }
    free(used); free(off); free(tried);
    return seq;
}

// Overwrites `count` random copies: a random source segment [s, s+len) is copied to a
// random target [t, t+len) with |s - t| >= len, so copies never overlap.
static void insert_repeats(Str *seq, int count, int len, Rng *rng) {
    const size_t n = seq->n;
    if (len < 0 || (size_t)len > n / 2) die("repeat length must be at most N/2");
    const size_t span = n - (size_t)len + 1; // valid start positions
    const size_t L = (size_t)len;
    for (int c = 0; c < count; ++c) {
        size_t s, below, above;
        do { // source s must leave a target: t in [0, s-L] or [s+L, span-1]
            s = rng_next(rng) % span;
            below = s >= L ? s - L + 1 : 0;
            above = s + L <= span - 1 ? span - (s + L) : 0;
        } while (below + above == 0);
        size_t r = rng_next(rng) % (below + above);
        size_t t = r < below ? r : s + L + (r - below);
        memmove(seq->s + t, seq->s + s, L);
    }
}

int main(int argc, char **argv) {
    size_t N = 10000;
    int k = 14, repeats = 0, repeat_len = 0;
    uint64_t seed = 1;
    const char *out_path = "data/ref_10000.txt";
    for (int a = 1; a + 1 < argc; a += 2) {
        const char *flag = argv[a], *val = argv[a + 1];
        if (!strcmp(flag, "--N")) N = strtoull(val, NULL, 10);
        else if (!strcmp(flag, "--k")) k = atoi(val);
        else if (!strcmp(flag, "--seed")) seed = strtoull(val, NULL, 10);
        else if (!strcmp(flag, "--repeats")) repeats = atoi(val);
        else if (!strcmp(flag, "--repeat-len")) repeat_len = atoi(val);
        else if (!strcmp(flag, "--out")) out_path = val;
        else { fprintf(stderr, "unknown flag %s\n", flag); return 1; }
    }
    Rng rng;
    rng_seed(&rng, seed);
    Str seq = gen_unique_kmer_seq(N, k, &rng);
    insert_repeats(&seq, repeats, repeat_len, &rng);
    FILE *f = fopen(out_path, "wb");
    if (!f) { fprintf(stderr, "cannot write %s\n", out_path); return 1; }
    fprintf(f, "%s\n", seq.s);
    fclose(f);
    return 0;
}
