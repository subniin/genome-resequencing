// SSAKE-style seed-and-extend assembler: sorted k-mer key indexes over the reads,
// extension in both directions.
// Build: gcc -O2 -std=c11 c/ssake.c -o build/c/ssake.exe
#include "common.h"

typedef struct { uint64_t key; uint32_t id; } KP;
typedef struct { uint64_t *keys; uint32_t *ids; size_t n; } KeyIndex;

static int kp_cmp(const void *pa, const void *pb) {
    const KP *a = pa, *b = pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return a->id < b->id ? -1 : a->id > b->id;
}

// Indexes the first kk bases (suffix == 0) or the last kk bases (suffix == 1) of each read.
static KeyIndex build_key_index(char **reads, size_t n, int L, int kk, int suffix) {
    const int end = suffix ? L - 1 : kk - 1; // index of the k-mer's last base
    KP *pairs = xrealloc(NULL, n * sizeof(KP));
    for (size_t i = 0; i < n; ++i) pairs[i] = (KP){kmer_code(reads[i], (size_t)end, kk), (uint32_t)i};
    qsort(pairs, n, sizeof(KP), kp_cmp);
    KeyIndex ix = {xrealloc(NULL, n * sizeof(uint64_t)), xrealloc(NULL, n * sizeof(uint32_t)), n};
    for (size_t i = 0; i < n; ++i) { ix.keys[i] = pairs[i].key; ix.ids[i] = pairs[i].id; }
    free(pairs);
    return ix;
}

// Block [lo, hi) of entries whose key equals `key`.
static void key_range(const KeyIndex *ix, uint64_t key, size_t *lo, size_t *hi) {
    size_t a = 0, b = ix->n;
    while (a < b) { size_t mid = (a + b) / 2; if (ix->keys[mid] < key) a = mid + 1; else b = mid; }
    *lo = a;
    b = ix->n;
    while (a < b) { size_t mid = (a + b) / 2; if (ix->keys[mid] <= key) a = mid + 1; else b = mid; }
    *hi = a;
}

// Finds the first unused read whose prefix (side 0) or suffix (side 1) overlaps the contig
// end (start) by k >= o bases with at most D mismatches, trying the largest k first.
// Returns the read id and sets *ov to k, or returns -1.
static long find_ext(const KeyIndex *ix, char **reads, const uint8_t *used, const Str *contig,
                     int side, int o, int kk, int L, int D, int *ov) {
    const int kmax = (int)contig->n < L - 1 ? (int)contig->n : L - 1;
    for (int k = kmax; k >= o; --k) {
        uint64_t key = side == 0 ? kmer_code(contig->s, contig->n - (size_t)k + (size_t)kk - 1, kk)
                                 : kmer_code(contig->s, (size_t)k - 1, kk);
        size_t lo, hi;
        key_range(ix, key, &lo, &hi);
        for (size_t x = lo; x < hi; ++x) {
            size_t id = ix->ids[x];
            if (used[id]) continue;
            int bad = side == 0 ? mismatches_bounded(reads[id], 0, contig->s, contig->n - (size_t)k, (size_t)k, D)
                                : mismatches_bounded(reads[id], (size_t)(L - k), contig->s, 0, (size_t)k, D);
            if (bad <= D) { *ov = k; return (long)id; }
        }
    }
    return -1;
}

// Assembles reads with seed-and-extend. Seeds are taken in input order; each seed is
// extended at its 3' end, then at its 5' end, using only unused reads. The first valid
// read at the largest overlap extends the contig. The index key is capped at 32 bases;
// the whole overlap is checked for mismatches.
static Str *assemble_ssake(char **reads, int n_reads, int L, int o, int D, size_t *n_out) {
    if (o < 1) die("min overlap must be >= 1");
    *n_out = 0;
    if (n_reads == 0) return NULL;
    if (o >= L) die("min overlap must be shorter than the reads");
    const size_t m = (size_t)dedupe_strings((const char **)reads, n_reads);
    const int kk = o < 32 ? o : 32;

    KeyIndex idx[2] = {build_key_index(reads, m, L, kk, 0), build_key_index(reads, m, L, kk, 1)};
    uint8_t *used = calloc(m, 1);
    Str *contigs = NULL;
    size_t nct = 0, capct = 0;

    for (size_t s = 0; s < m; ++s) {
        if (used[s]) continue;
        used[s] = 1;
        Str contig = {0};
        str_append(&contig, reads[s], (size_t)L);

        for (int side = 0; side < 2; ++side) { // 0: 3' end (prefix index), 1: 5' end (suffix index)
            int ov;
            long id;
            while ((id = find_ext(&idx[side], reads, used, &contig, side, o, kk, L, D, &ov)) >= 0) {
                used[id] = 1;
                if (side == 0) str_append(&contig, reads[id] + ov, (size_t)(L - ov));
                else str_prepend(&contig, reads[id], (size_t)(L - ov));
            }
        }
        GROW(contigs, nct, capct);
        contigs[nct++] = contig;
    }
    *n_out = nct;
    return contigs;
}

int main(int argc, char **argv) {
    const char *reads_path = NULL, *out_path = NULL;
    int o = 15, D = 2;
    for (int a = 1; a + 1 < argc; a += 2) {
        if (!strcmp(argv[a], "--reads")) reads_path = argv[a + 1];
        else if (!strcmp(argv[a], "--o")) o = atoi(argv[a + 1]);
        else if (!strcmp(argv[a], "--D")) D = atoi(argv[a + 1]);
        else if (!strcmp(argv[a], "--out")) out_path = argv[a + 1];
        else { fprintf(stderr, "unknown flag %s\n", argv[a]); return 1; }
    }
    if (!reads_path || !out_path) {
        fprintf(stderr, "usage: ssake --reads FILE --o O --D D --out FASTA\n");
        return 1;
    }
    // Reads: one per line, upper-cased, all of the first read's length, ACGT only.
    Lines ln = lines_load(reads_path);
    char **reads = xrealloc(NULL, (ln.n + 1) * sizeof(char *));
    int n = 0, L = 0;
    for (size_t i = 0; i < ln.n; ++i) {
        char *line = ln.v[i];
        size_t len = strlen(line);
        if (!len) continue;
        for (char *c = line; *c; ++c) {
            *c = (char)toupper((unsigned char)*c);
            if (*c != 'A' && *c != 'C' && *c != 'G' && *c != 'T') die("line %zu: non-ACGT base: %c", i + 1, *c);
        }
        if (n == 0) L = (int)len;
        else if ((int)len != L) die("line %zu: read length %zu differs from %d", i + 1, len, L);
        reads[n++] = line;
    }
    size_t nc;
    Str *contigs = assemble_ssake(reads, n, L, o, D, &nc);
    write_fasta(out_path, contigs, nc);
    return 0;
}
