// Benchmark: greedy overlap assembly with all-pairs overlap search.
// Build: gcc -O2 -std=c11 c/greedy_overlap.c -o build/c/greedy_overlap.exe
#include "common.h"

// Best overlap of a's suffix with b's prefix, length in [o, min(|a|,|b|)-1], at most D mismatches.
static int best_overlap(const char *a, size_t na, const char *b, size_t nb, int o, int D) {
    int kmax = (int)(na < nb ? na : nb) - 1;
    for (int k = kmax; k >= o; --k)
        if (mismatches_bounded(a, na - (size_t)k, b, 0, (size_t)k, D) <= D) return k;
    return 0;
}

static int uf_find(int *p, int x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
}

typedef struct { int k, i, j; } Cand;

static int cand_cmp(const void *pa, const void *pb) {
    const Cand *x = pa, *y = pb;
    if (x->k != y->k) return x->k > y->k ? -1 : 1;
    if (x->i != y->i) return x->i < y->i ? -1 : 1;
    return x->j < y->j ? -1 : x->j > y->j;
}

// Greedy assembly: all pairs scored once, links added from the longest overlap down. A link
// is kept when both ends are free and it does not close a cycle. Each chain is one contig.
static Str *assemble_greedy(const char **input, int n_in, int o, int D, size_t *n_contigs) {
    if (o < 1) die("min overlap must be >= 1");
    for (int i = 0; i < n_in; ++i)
        if (o >= (int)strlen(input[i])) die("min overlap must be shorter than the reads");

    const char **reads = input;
    const int m = dedupe_strings(reads, n_in);

    Cand *cand = NULL;
    size_t nc = 0, capc = 0;
    size_t *len = xrealloc(NULL, (size_t)m * sizeof(size_t));
    for (int i = 0; i < m; ++i) len[i] = strlen(reads[i]);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j) {
            if (i == j) continue;
            int k = best_overlap(reads[i], len[i], reads[j], len[j], o, D);
            if (k > 0) { GROW(cand, nc, capc); cand[nc++] = (Cand){k, i, j}; }
        }
    qsort(cand, nc, sizeof(Cand), cand_cmp);

    int *next = xrealloc(NULL, (size_t)m * sizeof(int)), *prev = xrealloc(NULL, (size_t)m * sizeof(int));
    int *ov = calloc((size_t)m + 1, sizeof(int)), *uf = xrealloc(NULL, (size_t)m * sizeof(int));
    for (int i = 0; i < m; ++i) { next[i] = prev[i] = -1; uf[i] = i; }
    for (size_t c = 0; c < nc; ++c) {
        int k = cand[c].k, i = cand[c].i, j = cand[c].j;
        if (next[i] != -1 || prev[j] != -1) continue;
        int ri = uf_find(uf, i), rj = uf_find(uf, j);
        if (ri == rj) continue; // would close a cycle
        next[i] = j; prev[j] = i; ov[i] = k; uf[ri] = rj;
    }

    Str *contigs = NULL;
    size_t cn = 0, ccap = 0;
    for (int s = 0; s < m; ++s) {
        if (prev[s] != -1) continue;
        Str t = {0};
        for (size_t x = 0; x < len[s]; ++x) str_push(&t, reads[s][x]);
        for (int cur = s; next[cur] != -1; cur = next[cur]) {
            int nx = next[cur];
            for (size_t x = (size_t)ov[cur]; x < len[nx]; ++x) str_push(&t, reads[nx][x]);
        }
        GROW(contigs, cn, ccap);
        contigs[cn++] = t;
    }
    free(cand); free(len); free(next); free(prev); free(ov); free(uf);
    *n_contigs = cn;
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
        fprintf(stderr, "usage: greedy_overlap --reads FILE --o O --D D --out FASTA\n");
        return 1;
    }
    Lines L = lines_load(reads_path);
    const char **reads = xrealloc(NULL, (L.n + 1) * sizeof(char *));
    int n = 0;
    for (size_t i = 0; i < L.n; ++i) if (L.v[i][0]) reads[n++] = L.v[i];
    size_t nc;
    Str *contigs = assemble_greedy(reads, n, o, D, &nc);
    write_fasta(out_path, contigs, nc);
    return 0;
}
