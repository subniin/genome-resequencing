// Generates a sample genome from a reference (SNPs only) and samples short reads from it.
// Build: gcc -O2 -std=c11 c/gen_reads.c -o build/c/gen_reads.exe
#include <math.h>
#include "common.h"
#ifdef _WIN32
#include <direct.h>
#define MKDIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MKDIR(p) mkdir(p, 0777)
#endif

typedef struct { size_t pos; char ref, alt; } Snp;

static const char kBase[4] = {'A', 'C', 'G', 'T'};

// Uniform integer in [0, urange], same method as libstdc++'s uniform_int_distribution.
static uint64_t uniform_int(Rng *r, uint64_t urange) {
    if (urange == UINT64_MAX) return rng_next(r);
    uint64_t range = urange + 1;
    unsigned __int128 prod = (unsigned __int128)rng_next(r) * range;
    uint64_t low = (uint64_t)prod;
    if (low < range) {
        uint64_t threshold = -range % range;
        while (low < threshold) {
            prod = (unsigned __int128)rng_next(r) * range;
            low = (uint64_t)prod;
        }
    }
    return (uint64_t)(prod >> 64);
}

// Uniform double in [0, 1), same as std::generate_canonical<double, 53>.
static double uniform_unit(Rng *r) {
    double v = (double)rng_next(r) / 18446744073709551616.0;
    return v >= 1.0 ? nextafter(1.0, 0.0) : v;
}

static void require_acgt(const Str *seq) {
    for (size_t i = 0; i < seq->n; ++i) {
        char c = seq->s[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') die("non-ACGT base in reference: %c", c);
    }
}

// Replaces `count` distinct positions of seq with a different base (positions drawn by rejection).
static Snp *insert_snps(Str *seq, int count, Rng *rng) {
    const size_t n = seq->n;
    uint8_t *chosen = calloc(n / 8 + 1, 1);
    Snp *snps = xrealloc(NULL, (size_t)count * sizeof(Snp));
    int got = 0;
    while (got < count) {
        size_t p = rng_next(rng) % n;
        if (chosen[p / 8] >> (p % 8) & 1) continue;
        chosen[p / 8] |= (uint8_t)(1 << (p % 8));
        char ref = seq->s[p], alt;
        do { alt = kBase[rng_next(rng) % 4]; } while (alt == ref);
        seq->s[p] = alt;
        snps[got++] = (Snp){p, ref, alt};
    }
    free(chosen);
    return snps;
}

static void mkdir_p(const char *path) {
    char *tmp = xrealloc(NULL, strlen(path) + 1);
    strcpy(tmp, path);
    for (char *c = tmp + 1; *c; ++c)
        if (*c == '/' || *c == '\\') { char s = *c; *c = 0; MKDIR(tmp); *c = s; }
    MKDIR(tmp);
    free(tmp);
}

int main(int argc, char **argv) {
    const char *ref_path = NULL, *out_dir = NULL;
    int L = 30, snps = 50;
    double coverage = 10, err = 0;
    uint64_t seed = 1;
    for (int a = 1; a + 1 < argc; a += 2) {
        const char *flag = argv[a], *val = argv[a + 1];
        if (!strcmp(flag, "--ref")) ref_path = val;
        else if (!strcmp(flag, "--L")) L = atoi(val);
        else if (!strcmp(flag, "--coverage")) coverage = atof(val);
        else if (!strcmp(flag, "--snps")) snps = atoi(val);
        else if (!strcmp(flag, "--err")) err = atof(val);
        else if (!strcmp(flag, "--seed")) seed = strtoull(val, NULL, 10);
        else if (!strcmp(flag, "--out")) out_dir = val;
        else { fprintf(stderr, "unknown flag %s\n", flag); return 1; }
    }
    if (!ref_path || !out_dir) {
        fprintf(stderr, "usage: gen_reads --ref FILE --L L --coverage C --snps K --err P --seed S --out DIR\n");
        return 1;
    }
    Rng rng;
    rng_seed(&rng, seed);
    Str sample = read_fasta_seq(ref_path);
    require_acgt(&sample);
    if (L <= 0 || (size_t)L > sample.n) die("read length must be in [1, N]");
    if (snps < 0 || (size_t)snps > sample.n) die("snp count must be in [0, N]");
    Snp *truth = insert_snps(&sample, snps, &rng);
    if (!(err >= 0 && err < 1)) die("err rate must be in [0, 1)"); // before any output
    mkdir_p(out_dir);

    char path[1024];
    snprintf(path, sizeof path, "%s/sample.txt", out_dir);
    FILE *f = fopen(path, "wb");
    if (!f) die("cannot write %s", path);
    fprintf(f, "%s\n", sample.s);
    fclose(f);

    snprintf(path, sizeof path, "%s/truth_snps.tsv", out_dir);
    f = fopen(path, "wb");
    if (!f) die("cannot write %s", path);
    for (int i = 0; i < snps; ++i) fprintf(f, "%zu\t%c\t%c\n", truth[i].pos, truth[i].ref, truth[i].alt);
    fclose(f);

    // M = round(N * coverage / L) reads at uniform starts in [0, N - L]; each base is
    // replaced by a different base with probability err.
    snprintf(path, sizeof path, "%s/reads.txt", out_dir);
    f = fopen(path, "wb");
    if (!f) die("cannot write %s", path);
    const size_t n = sample.n;
    const size_t m = (size_t)llround((double)n * coverage / L);
    char *read = xrealloc(NULL, (size_t)L + 1);
    for (size_t i = 0; i < m; ++i) {
        size_t st = (size_t)uniform_int(&rng, n - (size_t)L);
        memcpy(read, sample.s + st, (size_t)L);
        read[L] = 0;
        if (err > 0) {
            for (int j = 0; j < L; ++j) {
                if (uniform_unit(&rng) >= err) continue;
                int b = 0;
                while (kBase[b] != read[j]) ++b;
                read[j] = kBase[(b + 1 + rng_next(&rng) % 3) % 4];
            }
        }
        fprintf(f, "%s\n", read);
    }
    fclose(f);
    return 0;
}
