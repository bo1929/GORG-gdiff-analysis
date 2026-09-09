#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* Reconcile the two directions of a gdiff sample file into a single sample.
 *
 * Symmetric pairs (X,Y) with X>Y lexicographically act as the group key.
 * For each pair the "XY" direction rows (genome_a=X, genome_b=Y) and the
 * "YX" direction rows (genome_a=Y, genome_b=X) are two lists of window
 * samples. To reconcile:
 *   - sort each direction's rows by distance d ascending, dropping NaN d,
 *   - for rank i = 1..max(n_XY, n_YX), pick the row with the i-th smallest d
 *     from each direction when present; a number beats NaN; if both numbers,
 *     take the smaller (ties -> prefer XY).
 *   - the emitted row carries ALL columns of the chosen row, so the
 *     per-window fields (qid start end strand reference lr_bg lr_ub) of the
 *     selected direction are preserved (mirrors reconcile_directional).
 *
 * Output keeps the SAME 11-column schema as the input:
 *   config genome_a genome_b qid start end strand reference d lr_bg lr_ub
 * Rows per pair = max(#xy, #yx), i.e. from ~2n to ~n.
 *
 * Usage: reconcile-gdiff-sample <input.tsv> > <reconciled.tsv>
 */

typedef struct {
    char *qid, *start, *end, *strand, *reference;
    double d;
    char *lr_bg, *lr_ub;
} Row;

typedef struct { Row *v; size_t n, cap; } RowArr;

static void rowarr_push(RowArr *a, Row x) {
    if (a->n == a->cap) {
        a->cap = a->cap ? a->cap * 2 : 16;
        a->v = realloc(a->v, a->cap * sizeof(Row));
    }
    a->v[a->n++] = x;
}

typedef struct {
    uint64_t h;
    char *config, *X, *Y;
    RowArr xy, yx;
} Entry;

static Entry *tab;
static size_t tab_cap, tab_mask, tab_used;

static uint64_t fnv1a(const char *s) {
    uint64_t h = 1469598103934665603ULL;
    while (*s) { h ^= (unsigned char)*s++; h *= 1099511628211ULL; }
    return h;
}

static uint64_t mix3(uint64_t a, uint64_t b, uint64_t c) {
    uint64_t h = a;
    h ^= b + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h ^= c + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

static void tab_grow(void) {
    size_t newcap = tab_cap ? tab_cap * 2 : 1 << 16;
    Entry *nt = calloc(newcap, sizeof(Entry));
    size_t nmask = newcap - 1;
    for (size_t i = 0; i < tab_cap; i++) {
        if (!tab[i].config) continue;
        size_t j = tab[i].h & nmask;
        while (nt[j].config) j = (j + 1) & nmask;
        nt[j] = tab[i];
    }
    free(tab);
    tab = nt;
    tab_cap = newcap;
    tab_mask = nmask;
}

static Entry *get_entry(const char *cfg, const char *X, const char *Y, int create) {
    uint64_t h = mix3(fnv1a(cfg), fnv1a(X), fnv1a(Y));
    size_t i = h & tab_mask;
    while (tab[i].config) {
        if (tab[i].h == h && !strcmp(tab[i].config, cfg)
            && !strcmp(tab[i].X, X) && !strcmp(tab[i].Y, Y))
            return &tab[i];
        i = (i + 1) & tab_mask;
    }
    if (!create) return NULL;
    if ((tab_used + 1) * 4 >= tab_cap * 3) {
        tab_grow();
        i = h & tab_mask;
        while (tab[i].config) i = (i + 1) & tab_mask;
    }
    Entry *e = &tab[i];
    e->h = h;
    e->config = strdup(cfg);
    e->X = strdup(X);
    e->Y = strdup(Y);
    tab_used++;
    return e;
}

/* ascending by d; NaN sorts last */
static int cmp_row(const void *a, const void *b) {
    const Row *r1 = a, *r2 = b;
    int ax = isnan(r1->d), ay = isnan(r2->d);
    if (ax && ay) return 0;
    if (ax) return 1;
    if (ay) return -1;
    if (r1->d < r2->d) return -1;
    if (r1->d > r2->d) return 1;
    return 0;
}

static int cmp_entry(const void *a, const void *b) {
    const Entry *e1 = a, *e2 = b;
    int r;
    if ((r = strcmp(e1->config, e2->config))) return r;
    if ((r = strcmp(e1->X, e2->X))) return r;
    return strcmp(e1->Y, e2->Y);
}

static void sort_into(const RowArr *src, Row *dst) {
    for (size_t i = 0; i < src->n; i++) dst[i] = src->v[i];
    qsort(dst, src->n, sizeof(Row), cmp_row);
}

static void emit_row(const Entry *e, const Row *r) {
    if (isnan(r->d))
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tnan\t%s\t%s\n",
               e->config, e->X, e->Y, r->qid, r->start, r->end, r->strand,
               r->reference, r->lr_bg, r->lr_ub);
    else
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.6g\t%s\t%s\n",
               e->config, e->X, e->Y, r->qid, r->start, r->end, r->strand,
               r->reference, r->d, r->lr_bg, r->lr_ub);
}

static int split(char *line, char **f, int maxf) {
    int n = 0;
    char *p = line;
    while (n < maxf) {
        char *t = strchr(p, '\t');
        if (!t) { f[n++] = p; break; }
        *t = 0; f[n++] = p; p = t + 1;
    }
    return n;
}

static char *dup_or_null(const char *s) {
    return s ? strdup(s) : strdup(".");
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <gdiff-sample.tsv> > <reconciled.tsv>\n", argv[0]);
        return 1;
    }
    FILE *fh = fopen(argv[1], "r");
    if (!fh) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }

    tab_grow();
    char *line = NULL; size_t linelen = 0; ssize_t r;
    size_t lineno = 0;
    while ((r = getline(&line, &linelen, fh)) != -1) {
        lineno++;
        if (lineno == 1) continue;
        /* strip trailing newline so the last field has no "\n" embedded */
        if (r > 0 && line[r-1] == '\n') line[r-1] = '\0';
        char *f[11]; int nf = split(line, f, 11);
        if (nf < 9) continue;
        char *cfg = f[0], *ga = f[1], *gb = f[2];
        char *ds  = f[8];
        int c = strcmp(ga, gb);
        if (c == 0) continue;
        char *X = c > 0 ? ga : gb;
        char *Y = c > 0 ? gb : ga;
        double d;
        if (!*ds || !strcmp(ds, ".") || !strcmp(ds, "nan") || !strcmp(ds, "NA")) d = NAN;
        else {
            char *end;
            d = strtod(ds, &end);
            if (end == ds || *end) d = NAN;
        }
        Row row;
        row.qid = dup_or_null(nf > 3 ? f[3] : NULL);
        row.start = dup_or_null(nf > 4 ? f[4] : NULL);
        row.end = dup_or_null(nf > 5 ? f[5] : NULL);
        row.strand = dup_or_null(nf > 6 ? f[6] : NULL);
        row.reference = dup_or_null(nf > 7 ? f[7] : NULL);
        row.lr_bg = dup_or_null(nf > 9 ? f[9] : NULL);
        row.lr_ub = dup_or_null(nf > 10 ? f[10] : NULL);
        row.d = d;
        Entry *e = get_entry(cfg, X, Y, 1);
        rowarr_push(c > 0 ? &e->xy : &e->yx, row);
    }
    free(line);
    fclose(fh);

    printf("config\tgenome_a\tgenome_b\tqid\tstart\tend\tstrand\treference\td\tlr_bg\tlr_ub\n");

    Entry *out = malloc(tab_used * sizeof(Entry));
    size_t k = 0;
    for (size_t i = 0; i < tab_cap; i++)
        if (tab[i].config) out[k++] = tab[i];
    qsort(out, k, sizeof(Entry), cmp_entry);

    size_t n_out = 0;
    for (size_t i = 0; i < k; i++) {
        Entry *e = &out[i];
        size_t n1 = e->xy.n, n2 = e->yx.n, nmax = n1 > n2 ? n1 : n2;
        Row *s1 = malloc((n1 ? n1 : 1) * sizeof(Row));
        Row *s2 = malloc((n2 ? n2 : 1) * sizeof(Row));
        sort_into(&e->xy, s1);
        sort_into(&e->yx, s2);

        for (size_t j = 0; j < nmax; j++) {
            int has1 = j < n1 && !isnan(s1[j].d);
            int has2 = j < n2 && !isnan(s2[j].d);
            const Row *chosen = NULL;
            if (has1 && has2) chosen = (s1[j].d <= s2[j].d) ? &s1[j] : &s2[j];
            else if (has1) chosen = &s1[j];
            else if (has2) chosen = &s2[j];
            if (!chosen) {
                /* both NaN: emit an all-NA row */
                Row na = { strdup("."), strdup("."), strdup("."), strdup("."),
                           strdup("."), NAN, strdup("."), strdup(".") };
                emit_row(e, &na);
                n_out++;
                continue;
            }
            emit_row(e, chosen);
            n_out++;
        }
        free(s1); free(s2);
    }
    free(out);
    fprintf(stderr, "wrote %zu reconciled rows for %zu pairs\n", n_out, k);
    return 0;
}