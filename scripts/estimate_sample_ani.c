#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* estimate_sample_ani.c
 *
 * Reconciles the two directions of a gdiff sample file (same grouping as
 * reconcile-gdiff-sample.c), then for each symmetric pair (X,Y) computes a
 * pairwise distance estimate from the reconciled sample:
 *
 *   1. Reconciled sample: sorted per-direction rows merged rank-wise
 *      (min per rank, number beats NaN) -> up to max(#XY,#YX) rows.
 *   2. Filter: keep reconciled rows whose lr_ub >= --chi-sq (default 3.841);
 *      rows with lr_ub missing/NA are kept (no information to filter on).
 *   3. Decision:
 *        portion = (#non-NA rows kept by filter) / (#non-NA rows)
 *        if portion < --min-portion (default 0.66): report unfiltered mean
 *        else:                                   report filtered mean
 *   4. Output one row per pair with the distance plus extra statistics
 *      (columns after genome_b):
 *        distance  num_filtered  alternative_mean  num_NA
 *        max_unfiltered_distance  max_distance  num_lr_ub_zero
 *
 * Usage: estimate-sample-ani <gdiff-sample.tsv> [--chi-sq 3.841] [--min-portion 0.66]
 *   output: config genome_a genome_b distance num_filtered alternative_mean
 *           num_NA max_unfiltered_distance max_distance num_lr_ub_zero
 */

typedef struct {
    double d;     /* reconciled distance */
    double lr_ub; /* corresponding lr_ub (NaN if missing) */
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

static double bad2nan(const char *s) {
    if (!s || !*s) return NAN;
    if (!strcmp(s, ".") || !strcmp(s, "nan") || !strcmp(s, "NA") || !strcmp(s, "-")) return NAN;
    char *end;
    double v = strtod(s, &end);
    if (end == s || *end) return NAN;
    return v;
}

/* statistics for one pair */
typedef struct {
    double distance;
    long num_filtered;
    double alternative_mean;
    long num_NA;
    double max_unfiltered;
    double max_distance;
    long num_lr_ub_zero;
} PairStats;

static PairStats compute_stats(const Entry *e, double chi_sq, double min_portion) {
    PairStats st;
    st.distance = NAN; st.alternative_mean = NAN; st.max_unfiltered = NAN;
    st.max_distance = NAN; st.num_filtered = 0; st.num_NA = 0; st.num_lr_ub_zero = 0;

    size_t n1 = e->xy.n, n2 = e->yx.n, nmax = n1 > n2 ? n1 : n2;
    Row *s1 = malloc((n1 ? n1 : 1) * sizeof(Row));
    Row *s2 = malloc((n2 ? n2 : 1) * sizeof(Row));
    sort_into(&e->xy, s1);
    sort_into(&e->yx, s2);

    /* reconcile: min per rank, number beats NaN */
    size_t i;
    double unf_mean_sum = 0, fil_mean_sum = 0;
    size_t unf_cnt = 0, fil_cnt = 0;
    for (i = 0; i < nmax; i++) {
        int has1 = i < n1 && !isnan(s1[i].d);
        int has2 = i < n2 && !isnan(s2[i].d);
        Row chosen;
        if (has1 && has2) chosen = (s1[i].d <= s2[i].d) ? s1[i] : s2[i];
        else if (has1) chosen = s1[i];
        else if (has2) chosen = s2[i];
        else { st.num_NA++; continue; }   /* reconciled row is NA */

        if (isnan(chosen.d)) { st.num_NA++; continue; }

        if (chosen.lr_ub == 0.0) st.num_lr_ub_zero++;

        /* unfiltered accumulator */
        unf_mean_sum += chosen.d; unf_cnt++;
        if (isnan(st.max_unfiltered) || chosen.d > st.max_unfiltered)
            st.max_unfiltered = chosen.d;

        /* filtered: keep if lr_ub >= chi_sq (or lr_ub missing -> keep) */
        int keep = isnan(chosen.lr_ub) || chosen.lr_ub >= chi_sq;
        if (keep) {
            fil_mean_sum += chosen.d; fil_cnt++;
            if (isnan(st.max_distance) || chosen.d > st.max_distance)
                st.max_distance = chosen.d;
        } else {
            st.num_filtered++;
        }
    }
    free(s1); free(s2);

    size_t non_na = unf_cnt;
    double portion = non_na > 0 ? (double)fil_cnt / (double)non_na : 0.0;
    double unf_mean = unf_cnt ? unf_mean_sum / unf_cnt : NAN;
    double fil_mean = fil_cnt ? fil_mean_sum / fil_cnt : NAN;

    if (portion < min_portion) {
        st.distance = unf_mean;
        st.alternative_mean = fil_mean;
    } else {
        st.distance = fil_mean;
        st.alternative_mean = unf_mean;
    }
    return st;
}

static void emit_stats(const Entry *e, const PairStats *st) {
    printf("%s\t%s\t%s", e->config, e->X, e->Y);
    printf(isnan(st->distance) ? "\tnan" : "\t%.6g", st->distance);
    printf("\t%ld", st->num_filtered);
    printf(isnan(st->alternative_mean) ? "\tnan" : "\t%.6g", st->alternative_mean);
    printf("\t%ld", st->num_NA);
    printf(isnan(st->max_unfiltered) ? "\tnan" : "\t%.6g", st->max_unfiltered);
    printf(isnan(st->max_distance) ? "\tnan" : "\t%.6g", st->max_distance);
    printf("\t%ld\n", st->num_lr_ub_zero);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <gdiff-sample.tsv> [--chi-sq 3.841] [--min-portion 0.66]\n", argv[0]);
        return 1;
    }
    const char *path = argv[1];
    double chi_sq = 3.841, min_portion = 0.66;
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--chi-sq") == 0 && i + 1 < argc) chi_sq = atof(argv[++i]);
        else if (strcmp(argv[i], "--min-portion") == 0 && i + 1 < argc) min_portion = atof(argv[++i]);
        else { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 1; }
    }

    FILE *fh = fopen(path, "r");
    if (!fh) { fprintf(stderr, "cannot open %s\n", path); return 1; }

    tab_grow();
    char *line = NULL; size_t linelen = 0; ssize_t nrl;
    size_t lineno = 0;
    while ((nrl = getline(&line, &linelen, fh)) != -1) {
        lineno++;
        if (lineno == 1) continue;
        if (nrl > 0 && line[nrl-1] == '\n') line[nrl-1] = '\0';
        char *f[11]; int nf = split(line, f, 11);
        if (nf < 9) continue;
        char *cfg = f[0], *ga = f[1], *gb = f[2];
        char *ds  = f[8];
        int c = strcmp(ga, gb);
        if (c == 0) continue;
        char *X = c > 0 ? ga : gb;
        char *Y = c > 0 ? gb : ga;
        Row row;
        row.d = bad2nan(ds);
        row.lr_ub = bad2nan(nf > 10 ? f[10] : NULL);
        Entry *e = get_entry(cfg, X, Y, 1);
        rowarr_push(c > 0 ? &e->xy : &e->yx, row);
    }
    free(line);
    fclose(fh);

    printf("config\tgenome_a\tgenome_b\tdistance\tnum_filtered\talternative_mean\tnum_NA\tmax_unfiltered_distance\tmax_distance\tnum_lr_ub_zero\n");

    Entry *out = malloc(tab_used * sizeof(Entry));
    size_t k = 0;
    for (size_t i = 0; i < tab_cap; i++)
        if (tab[i].config) out[k++] = tab[i];
    qsort(out, k, sizeof(Entry), cmp_entry);

    for (size_t i = 0; i < k; i++) {
        PairStats st = compute_stats(&out[i], chi_sq, min_portion);
        emit_stats(&out[i], &st);
    }
    free(out);
    fprintf(stderr, "wrote %zu pairs\n", k);
    return 0;
}