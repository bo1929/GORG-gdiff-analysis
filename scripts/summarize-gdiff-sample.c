#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

typedef struct { double *v; size_t n, cap; } Arr;

static void arr_push(Arr *a, double x) {
    if (a->n == a->cap) {
        a->cap = a->cap ? a->cap * 2 : 16;
        a->v = realloc(a->v, a->cap * sizeof(double));
    }
    a->v[a->n++] = x;
}

typedef struct {
    uint64_t h;
    char *config, *X, *Y;
    Arr xy, yx;
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

static int cmp_nan_last(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    int ax = isnan(x), ay = isnan(y);
    if (ax && ay) return 0;
    if (ax) return 1;
    if (ay) return -1;
    return (x > y) - (x < y);
}

static double mean_of(const Arr *a, size_t *nnan) {
    double s = 0; size_t k = 0, c = 0;
    for (size_t i = 0; i < a->n; i++)
        if (isnan(a->v[i])) c++;
        else { s += a->v[i]; k++; }
    if (nnan) *nnan = c;
    return k ? s / k : NAN;
}

static double median_of(const Arr *a, size_t *nnan) {
    size_t k = 0, c = 0;
    double *tmp = malloc((a->n ? a->n : 1) * sizeof(double));
    for (size_t i = 0; i < a->n; i++)
        if (isnan(a->v[i])) c++;
        else tmp[k++] = a->v[i];
    if (nnan) *nnan = c;
    double m = NAN;
    if (k) {
        qsort(tmp, k, sizeof(double), cmp_nan_last);
        m = (k & 1) ? tmp[k / 2] : (tmp[k / 2 - 1] + tmp[k / 2]) / 2.0;
    }
    free(tmp);
    return m;
}

static double reconcile_mean(const Arr *xy, const Arr *yx, size_t *nnan) {
    size_t n1 = xy->n, n2 = yx->n, nmax = n1 > n2 ? n1 : n2;
    if (!nmax) nmax = 1;
    double *s1 = malloc((n1 ? n1 : 1) * sizeof(double));
    double *s2 = malloc((n2 ? n2 : 1) * sizeof(double));
    memcpy(s1, xy->v, n1 * sizeof(double));
    memcpy(s2, yx->v, n2 * sizeof(double));
    qsort(s1, n1, sizeof(double), cmp_nan_last);
    qsort(s2, n2, sizeof(double), cmp_nan_last);
    double sum = 0; size_t cnt = 0, c = 0;
    for (size_t i = 0; i < nmax; i++) {
        double v1 = i < n1 ? s1[i] : NAN;
        double v2 = i < n2 ? s2[i] : NAN;
        double ch;
        if (isnan(v1) && isnan(v2)) ch = NAN;
        else if (isnan(v1)) ch = v2;
        else if (isnan(v2)) ch = v1;
        else ch = v1 < v2 ? v1 : v2;
        if (isnan(ch)) c++;
        else { sum += ch; cnt++; }
    }
    free(s1); free(s2);
    if (nnan) *nnan = c;
    return cnt ? sum / cnt : NAN;
}

static void emit(double v, size_t nan) {
    if (isnan(v)) printf("\tnan\t%zu", nan);
    else printf("\t%.6f\t%zu", v, nan);
}

static int cmp_entry(const void *a, const void *b) {
    const Entry *e1 = a, *e2 = b;
    int r;
    if ((r = strcmp(e1->config, e2->config))) return r;
    if ((r = strcmp(e1->X, e2->X))) return r;
    return strcmp(e1->Y, e2->Y);
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

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <gdiff-sample.tsv>\n", argv[0]);
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
        char *f[11]; int nf = split(line, f, 11);
        if (nf < 9) continue;
        char *cfg = f[0], *ga = f[1], *gb = f[2], *ds = f[8];
        int c = strcmp(ga, gb);
        if (c == 0) continue;
        char *X = c > 0 ? ga : gb;
        char *Y = c > 0 ? gb : ga;
        double d;
        if (!*ds || !strcmp(ds, ".")) d = NAN;
        else {
            char *end;
            d = strtod(ds, &end);
            if (end == ds || *end) d = NAN;
        }
        Entry *e = get_entry(cfg, X, Y, 1);
        arr_push(c > 0 ? &e->xy : &e->yx, d);
    }
    free(line);
    fclose(fh);

    printf("config\tgenome_a\tgenome_b"
           "\tx\tna_x"
           "\ty\tna_y"
           "\tz\tna_z\n");

    /* collect used entries and sort deterministically */
    Entry *out = malloc(tab_used * sizeof(Entry));
    size_t k = 0;
    for (size_t i = 0; i < tab_cap; i++)
        if (tab[i].config) out[k++] = tab[i];
    qsort(out, k, sizeof(Entry), cmp_entry);

    for (size_t i = 0; i < k; i++) {
        Entry *e = &out[i];
        size_t n1, n2;
        double m1 = mean_of(&e->xy, &n1);
        double m2 = mean_of(&e->yx, &n2);
        size_t na_min = n1 + n2;
        double symmin = (isnan(m1) || isnan(m2)) ? NAN : (m1 < m2 ? m1 : m2);

        double md1 = median_of(&e->xy, &n1);
        double md2 = median_of(&e->yx, &n2);
        size_t na_med = n1 + n2;
        double symmed = (isnan(md1) || isnan(md2)) ? NAN : (md1 + md2) / 2.0;

        size_t na_rc;
        double rc = reconcile_mean(&e->xy, &e->yx, &na_rc);

        printf("%s\t%s\t%s", e->config, e->X, e->Y);
        emit(symmin, na_min);
        emit(symmed, na_med);
        emit(rc, na_rc);
        printf("\n");
    }
    free(out);
    return 0;
}
