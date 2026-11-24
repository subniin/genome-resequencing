// Shared helpers: errors, strings, file reading, RNG, read dedupe.
#ifndef COMMON_H
#define COMMON_H
#define __USE_MINGW_ANSI_STDIO 1
#include <ctype.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void die(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fputs("error: ", stderr);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    exit(1);
}

static void *xrealloc(void *p, size_t n) {
    p = realloc(p, n ? n : 1);
    if (!p) die("out of memory");
    return p;
}

// Makes room for element index n in array p (capacity cap).
#define GROW(p, n, cap) do { if ((n) >= (cap)) { (cap) = (cap) ? (cap) * 2 : 16; \
    (p) = xrealloc((p), (cap) * sizeof(*(p))); } } while (0)

// Growable NUL-terminated string. A zeroed Str is empty (s may be NULL until first push).
typedef struct { char *s; size_t n, cap; } Str;

static void str_push(Str *t, char c) {
    if (t->n + 1 >= t->cap) {
        t->cap = t->cap ? t->cap * 2 : 64;
        t->s = xrealloc(t->s, t->cap);
    }
    t->s[t->n++] = c;
    t->s[t->n] = 0;
}

static void str_append(Str *t, const char *p, size_t len) {
    for (size_t i = 0; i < len; ++i) str_push(t, p[i]);
}

// Inserts p[0..len) in front of t.
static void str_prepend(Str *t, const char *p, size_t len) {
    size_t old = t->n;
    str_append(t, p, len); // grows t
    memmove(t->s + len, t->s, old);
    memcpy(t->s, p, len);
}

static char *read_file(const char *path, size_t *len) {
    FILE *f = fopen(path, "rb");
    if (!f) die("cannot open %s", path);
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buf = xrealloc(NULL, (size_t)sz + 1);
    size_t got = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    buf[got] = 0;
    if (len) *len = got;
    return buf;
}

// All lines of a file (CR/LF stripped, empty lines kept so index + 1 is the line number).
typedef struct { char *buf; char **v; size_t n; } Lines;

static Lines lines_load(const char *path) {
    Lines L = {0};
    size_t len, cap = 0;
    L.buf = read_file(path, &len);
    char *p = L.buf, *end = L.buf + len;
    while (p < end) {
        char *nl = memchr(p, '\n', (size_t)(end - p));
        char *e = nl ? nl : end;
        char *next = nl ? nl + 1 : end;
        if (e > p && e[-1] == '\r') e--;
        *e = 0;
        GROW(L.v, L.n, cap);
        L.v[L.n++] = p;
        p = next;
    }
    return L;
}

// Sequence of a FASTA (or plain text) file: skips '>' lines, upper-cases bases.
static Str read_fasta_seq(const char *path) {
    Lines L = lines_load(path);
    Str seq = {0};
    for (size_t i = 0; i < L.n; ++i) {
        if (L.v[i][0] == '>') continue;
        for (char *c = L.v[i]; *c; ++c) str_push(&seq, (char)toupper((unsigned char)*c));
    }
    if (!seq.s) str_push(&seq, 'A'), seq.n = 0, seq.s[0] = 0;
    free(L.v); free(L.buf);
    return seq;
}

static void write_fasta(const char *path, const Str *contigs, size_t n) {
    FILE *f = fopen(path, "wb");
    if (!f) die("cannot write %s", path);
    for (size_t i = 0; i < n; ++i) fprintf(f, ">contig_%zu\n%s\n", i, contigs[i].s);
    fclose(f);
}

// Counts mismatches in a[ai..ai+len) vs b[bi..bi+len); stops at limit+1 once above limit.
static int mismatches_bounded(const char *a, size_t ai, const char *b, size_t bi, size_t len, int limit) {
    int mm = 0;
    for (size_t i = 0; i < len; ++i)
        if (a[ai + i] != b[bi + i] && ++mm > limit) return limit + 1;
    return mm;
}

// mt19937_64, same stream as std::mt19937_64.
typedef struct { uint64_t mt[312]; int idx; } Rng;

static void rng_seed(Rng *r, uint64_t seed) {
    r->mt[0] = seed;
    for (int i = 1; i < 312; ++i)
        r->mt[i] = 6364136223846793005ULL * (r->mt[i - 1] ^ (r->mt[i - 1] >> 62)) + (uint64_t)i;
    r->idx = 312;
}

static uint64_t rng_next(Rng *r) {
    if (r->idx >= 312) {
        for (int i = 0; i < 312; ++i) {
            uint64_t x = (r->mt[i] & 0xFFFFFFFF80000000ULL) | (r->mt[(i + 1) % 312] & 0x7FFFFFFFULL);
            uint64_t xa = x >> 1;
            if (x & 1) xa ^= 0xB5026F5AA96619E9ULL;
            r->mt[i] = r->mt[(i + 156) % 312] ^ xa;
        }
        r->idx = 0;
    }
    uint64_t y = r->mt[r->idx++];
    y ^= (y >> 29) & 0x5555555555555555ULL;
    y ^= (y << 17) & 0x71D67FFFEDA60000ULL;
    y ^= (y << 37) & 0xFFF7EEE000000000ULL;
    y ^= y >> 43;
    return y;
}

static int base_index(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        default:  return 3; // 'T'
    }
}

// Code of the k-mer that ends at index `end` (end + 1 >= k, k in [1, 32]).
static uint64_t kmer_code(const char *seq, size_t end, int k) {
    uint64_t code = 0;
    for (size_t t = end + 1 - (size_t)k; t <= end; ++t) code = (code << 2) | (uint64_t)base_index(seq[t]);
    return code;
}

static const char **g_strs; // comparator context for qsort

static int str_idx_cmp(const void *pa, const void *pb) {
    int a = *(const int *)pa, b = *(const int *)pb;
    int c = strcmp(g_strs[a], g_strs[b]);
    return c ? c : (a < b ? -1 : a > b);
}

// Keeps the first copy of each distinct string, in input order. Returns the new count.
static int dedupe_strings(const char **v, int n) {
    int *ord = xrealloc(NULL, (size_t)n * sizeof(int));
    uint8_t *keep = calloc((size_t)n + 1, 1);
    for (int i = 0; i < n; ++i) ord[i] = i;
    g_strs = v;
    qsort(ord, (size_t)n, sizeof(int), str_idx_cmp);
    for (int i = 0; i < n; ++i)
        if (i == 0 || strcmp(v[ord[i]], v[ord[i - 1]]) != 0) keep[ord[i]] = 1;
    int m = 0;
    for (int i = 0; i < n; ++i) if (keep[i]) v[m++] = v[i];
    free(ord); free(keep);
    return m;
}
#endif
