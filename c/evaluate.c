// Evaluates contigs against the reference (SNP calls) and the sample (accuracy).
// Build: gcc -O2 -std=c11 c/evaluate.c -o build/c/evaluate.exe
#include "common.h"

typedef struct { uint64_t code, pos; } RefKmer;
typedef struct { int cpos; int64_t rpos; } Anchor;

static int refk_cmp(const void *pa, const void *pb) {
    const RefKmer *a = pa, *b = pb;
    if (a->code != b->code) return a->code < b->code ? -1 : 1;
    return a->pos < b->pos ? -1 : a->pos > b->pos;
}

static int size_desc(const void *pa, const void *pb) {
    size_t a = *(const size_t *)pa, b = *(const size_t *)pb;
    return a > b ? -1 : a < b;
}

static int64_t diag(const Anchor *x) { return x->rpos - x->cpos; }

int main(int argc, char **argv) {
    const char *contigs_path = NULL, *ref_path = NULL, *sample_path = NULL, *truth_path = NULL;
    int k = 24;
    for (int a = 1; a + 1 < argc; a += 2) {
        const char *flag = argv[a], *val = argv[a + 1];
        if (!strcmp(flag, "--contigs")) contigs_path = val;
        else if (!strcmp(flag, "--ref")) ref_path = val;
        else if (!strcmp(flag, "--sample")) sample_path = val;
        else if (!strcmp(flag, "--truth")) truth_path = val;
        else if (!strcmp(flag, "--k")) k = atoi(val);
        else { fprintf(stderr, "unknown flag %s\n", flag); return 1; }
    }
    if (!contigs_path || !ref_path || !sample_path || !truth_path) {
        fprintf(stderr, "usage: evaluate --contigs F --ref F --sample F --truth TSV [--k K (default 24)]\n");
        return 1;
    }
    if (k < 1 || k > 32) die("k must be in [1, 32]");

    // Contigs: FASTA records, sequence lines joined.
    Lines cl = lines_load(contigs_path);
    Str *contigs = NULL;
    size_t nc = 0, capc = 0;
    for (size_t i = 0; i < cl.n; ++i) {
        if (cl.v[i][0] == '>') {
            GROW(contigs, nc, capc);
            memset(&contigs[nc++], 0, sizeof(Str));
        } else if (nc) {
            for (char *c = cl.v[i]; *c; ++c) str_push(&contigs[nc - 1], *c);
        }
    }
    for (size_t i = 0; i < nc; ++i) if (!contigs[i].s) str_push(&contigs[i], 'A'), contigs[i].n = 0, contigs[i].s[0] = 0;

    Str ref = read_fasta_seq(ref_path), sample = read_fasta_seq(sample_path);
    const size_t N = ref.n;

    // Truth SNPs: truth_alt[pos] is the planted base (0 = none).
    char *truth_alt = calloc(N + 1, 1);
    size_t n_truth = 0;
    Lines tl = lines_load(truth_path);
    for (size_t i = 0; i < tl.n; ++i) {
        char *line = tl.v[i];
        if (!line[0]) continue;
        size_t p = strtoull(line, NULL, 10);
        char *t2 = strrchr(line, '\t');
        if (p < N) truth_alt[p] = t2[1];
        n_truth++;
    }

    const uint64_t mask = k == 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);
    // Unique reference k-mers as (code, position), sorted by code; repeated codes are dropped.
    RefKmer *refK = xrealloc(NULL, (N >= (size_t)k ? N - (size_t)k + 1 : 0) * sizeof(RefKmer));
    size_t nk = 0;
    uint64_t code = 0;
    for (size_t q = 0; q < N; ++q) {
        code = ((code << 2) | (uint64_t)base_index(ref.s[q])) & mask;
        if (q + 1 >= (size_t)k) refK[nk++] = (RefKmer){code, q + 1 - (size_t)k};
    }
    qsort(refK, nk, sizeof(RefKmer), refk_cmp);
    size_t kept = 0;
    for (size_t i = 0; i < nk;) {
        size_t j = i + 1;
        while (j < nk && refK[j].code == refK[i].code) ++j;
        if (j - i == 1) refK[kept++] = refK[i];
        i = j;
    }
    nk = kept;

    size_t *lens = xrealloc(NULL, (nc + 1) * sizeof(size_t));
    size_t totalBases = 0, alignedBases = 0, unalignedBases = 0, matchedBases = 0;
    uint8_t *fl = calloc(N + 1, 1); // bit 0: covered, bit 1: SNP already counted
    size_t called = 0, hits = 0;
    int misassemblies = 0, unaligned_contigs = 0;
    Anchor *anchors = NULL;
    size_t capa = 0;

    for (size_t ci = 0; ci < nc; ++ci) {
        const char *c = contigs[ci].s;
        const int cn = (int)contigs[ci].n;
        lens[ci] = (size_t)cn;
        totalBases += (size_t)cn;

        size_t na = 0;
        uint64_t cc = 0;
        for (int i = 0; i < cn; ++i) {
            cc = ((cc << 2) | (uint64_t)base_index(c[i])) & mask;
            if (i + 1 < k) continue;
            size_t lo = 0, hi = nk; // lower_bound on code
            while (lo < hi) { size_t mid = (lo + hi) / 2; if (refK[mid].code < cc) lo = mid + 1; else hi = mid; }
            if (lo < nk && refK[lo].code == cc) {
                GROW(anchors, na, capa);
                anchors[na++] = (Anchor){i + 1 - k, (int64_t)refK[lo].pos};
            }
        }
        // Drop isolated spurious anchors: diagonal differs from both neighbours.
        if (na > 1) {
            size_t w = 0;
            Anchor *tmp = xrealloc(NULL, na * sizeof(Anchor));
            for (size_t i = 0; i < na; ++i) {
                int prevOk = i > 0 && diag(&anchors[i]) == diag(&anchors[i - 1]);
                int nextOk = i + 1 < na && diag(&anchors[i]) == diag(&anchors[i + 1]);
                if (prevOk || nextOk) tmp[w++] = anchors[i];
            }
            memcpy(anchors, tmp, w * sizeof(Anchor));
            free(tmp);
            na = w;
        }
        if (na == 0) { unaligned_contigs++; unalignedBases += (size_t)cn; continue; }

        alignedBases += (size_t)cn;
        // Split anchors into segments of constant diagonal; a contig base is scored once.
        size_t start = 0;
        int prevEnd = 0;
        for (size_t a = 1; a <= na; ++a) {
            int split = a == na || diag(&anchors[a]) != diag(&anchors[a - 1]);
            if (!split) continue;
            if (a < na) misassemblies++;

            int64_t d = diag(&anchors[start]);
            int cFrom = start == 0 ? 0 : anchors[start].cpos;
            if (prevEnd > cFrom) cFrom = prevEnd;
            int cTo = a == na ? cn : anchors[a - 1].cpos + k; // exclusive
            prevEnd = cTo;
            for (int i = cFrom; i < cTo; ++i) {
                int64_t q = i + d;
                if (q < 0 || q >= (int64_t)N) continue;
                fl[q] |= 1;
                if (c[i] == sample.s[q]) matchedBases++;
                if (c[i] != ref.s[q] && !(fl[q] & 2)) { // called SNP at reference position q
                    fl[q] |= 2;
                    called++;
                    if (truth_alt[q] == c[i]) hits++;
                }
            }
            start = a;
        }
    }

    size_t n50 = 0, total = 0, acc = 0, covered = 0;
    qsort(lens, nc, sizeof(size_t), size_desc);
    for (size_t i = 0; i < nc; ++i) total += lens[i];
    for (size_t i = 0; i < nc; ++i) {
        acc += lens[i];
        if (acc * 2 >= total) { n50 = lens[i]; break; }
    }
    for (size_t i = 0; i < N; ++i) covered += fl[i] & 1;

    printf("accuracy=%.6f\nsnp_recall=%.6f\nsnp_precision=%.6f\nn_contigs=%d\nn50=%zu\ncovered_frac=%.6f\n"
           "misassemblies=%d\nunaligned_contigs=%d\nunaligned_bases_frac=%.6f\n",
           alignedBases ? (double)matchedBases / alignedBases : 0.0,
           n_truth ? (double)hits / n_truth : 0.0,
           called ? (double)hits / called : 0.0,
           (int)nc, n50, N ? (double)covered / N : 0.0, misassemblies, unaligned_contigs,
           totalBases ? (double)unalignedBases / totalBases : 0.0);
    return 0;
}
