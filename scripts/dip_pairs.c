/* dip_pairs.c -- Hartigan's dip statistic per symmetric pair of a gdiff sample file.
 *
 * Pipeline (identical to scripts/estimate_sample_ani.c for the sample handling):
 *
 *   1. Symmetric key: (config, X, Y) with X > Y lexicographically.  Rows with
 *      genome_a == genome_b are skipped.  Rows where genome_a < genome_b go into
 *      the "yx" direction, the others into "xy".
 *   2. Reconcile the two directions: sort each direction's rows by distance d
 *      ascending (NaN sorts last), then for rank i = 1..max(n_xy, n_yx) take the
 *      i-th smallest d of each direction; a number beats NaN, and if both are
 *      numbers the smaller one wins (ties -> xy).  The reconciled array is
 *      therefore already sorted ascending and needs no second sort.
 *   3. Filter: a reconciled row is kept iff its lr_ub is present and
 *      lr_ub > chi_sq (--chi-sq, default 3.841).  Missing/NA lr_ub rows are
 *      dropped, exactly as in estimate_sample_ani.c.
 *   4. Selection: portion = #kept / #non-NA reconciled rows.
 *        portion >  min_portion  -> work on the FILTERED set  (--min-portion, 0.66)
 *        portion <= min_portion  -> work on the UNFILTERED set
 *      This mirrors the iddd switch of estimate_sample_ani.c and select_sample()
 *      in simulations/cv_analysis.R.  The unfiltered set is the full reconciled
 *      sample of finite d.
 *   5. dip = Hartigan & Hartigan dip statistic of the selected set, computed with
 *      the greatest-convex-minorant / least-concave-majorant algorithm of
 *      ALGORITHM AS 217, Appl. Statist. (1985) 34(3):320-325 -- the same
 *      algorithm as R's diptest::dip().  This file is a dependency-free C port of
 *      diptest's src/dip.c (P. M. Hartigan, f2c'd; fixes and speedups by
 *      Martin Maechler), validated against diptest 0.77-2 (see validate_dip.R).
 *
 * Output (TSV on stdout, one row per symmetric pair, sorted by
 * config,genome_a,genome_b):
 *
 *   config genome_a genome_b dip n_dip portion n_kept n_unfiltered used
 *   num_filtered num_NA max_unfiltered_distance max_dip_distance num_lr_ub_zero
 *
 *   dip, n_dip            the statistic and the size of the set it came from
 *   portion               n_kept / n_unfiltered
 *   used                  "filtered" | "unfiltered" (which set dip came from)
 *   max_dip_distance      max d over the selected set
 *   the remaining columns identical in meaning to estimate-sample-ani's output
 *
 * Usage:
 *   dip-pairs <gdiff-sample.tsv> [--chi-sq 3.841] [--min-portion 0.66]
 *                                [--min-n 4] [--jobs N] [--stream|--no-stream]
 *
 *   --min-n N        pairs whose selected set is smaller than N points (or < 2)
 *                    get dip = nan.  Default 2, matching diptest::dip(), which
 *                    is defined for n >= 2; raise it to ignore tiny samples.
 *   --jobs N         threads; 0 = all cores, 1 = serial (default 0)
 *   --stream         assume the rows of each pair are contiguous and emit as we
 *                    go (constant memory).  Auto-enabled for inputs >= 256 MiB.
 *                    Rows are emitted in input order, so the pair order follows
 *                    the file; content is identical to --no-stream.
 *   --no-stream      always buffer everything, then emit sorted by
 *                    (config, genome_a, genome_b).  Peak memory is ~16 bytes
 *                    per buffered input row, so prefer --stream for multi-GB
 *                    files (a 1.8 GB / 250M-row sample needs ~4 GB that way).
 *   --dump-selected F  also write "config<TAB>X<TAB>Y<TAB>d,d,d,..." for the
 *                    selected set of every pair to F (validation/debugging)
 *   --precision N    significant digits for floating point output (default 10)
 *
 * Build:  make -C scripts                (serial if no OpenMP)
 *         cc -O3 -fopenmp -o dip-pairs scripts/dip_pairs.c -lm
 * macOS + Homebrew libomp:
 *         cc -O3 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
 *            -L/opt/homebrew/opt/libomp/lib -lomp -o dip-pairs \
 *            scripts/dip_pairs.c -lm
 */

/*
 * We use getline(), ssize_t and off_t, which ISO C does not declare.  Ask for
 * POSIX.1-2008 plus the default (BSD/glibc) namespace *before* any header is
 * included; otherwise a strict -std=c99 build on x86-64 Linux/glibc fails with
 * "implicit declaration of getline" and "unknown type name ssize_t/off_t"
 * (macOS headers happen to expose all of these, so this only shows up there).
 * _DEFAULT_SOURCE is unnecessary on glibc >= 2.12 but is harmless and covers
 * older glibc, where strdup() sat outside _POSIX_C_SOURCE.
 *
 * strdup() itself is deliberately not used: Darwin gates it on _DARWIN_C_SOURCE
 * even when _POSIX_C_SOURCE is set, so xstrdup() below is open-coded instead.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE 1
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <limits.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define STREAM_THRESHOLD ((off_t)268435456)   /* 256 MiB */

/* ================================================================== */
/* Hartigan & Hartigan dip statistic (ALGORITHM AS 217)               */
/* ================================================================== */

/* Workspace for dip(): 5*(n+1) ints.  Callers normally let dip() allocate. */
typedef struct {
    int *mn, *mj, *gcm, *lcm;
    size_t n;
} DipWorkspace;

static void dip_ws_free(DipWorkspace *w)
{
    free(w->mn); free(w->mj); free(w->gcm); free(w->lcm);
    w->mn = w->mj = w->gcm = w->lcm = NULL;
    w->n = 0;
}

/* Hartigan & Hartigan dip statistic of a pre-sorted array.
 *
 *   x   ascending sample (ties allowed), length n
 *   ws  workspace; pass NULL to allocate/free internally
 *   lo/hi  optional 1-based modal-interval indices
 *   min_is_0  0 -> a flat sample yields the theoretical floor 2/(2n);
 *             1 -> it yields 0 (diptest's `min.is.0`)
 *
 * Returns NAN for n < 1 or for an unsorted x.  O(n) time for unimodal input,
 * and the value is exactly diptest::dip() to the last ulp for n <= 2^31.
 */
static double dip(const double *x, int n, DipWorkspace *ws, int *lo, int *hi,
                  int min_is_0)
{
    int low, high, l_gcm, l_lcm;
    int mnj, mnmnj, mjk, mjmjk, ig = 1, ih = 1, iv, ix, i, j, k;
    double dip_l, dip_u, dipnew, dip_val;
    DipWorkspace w, *W;

    if (n < 1) return NAN;
    for (k = 2; k <= n; ++k)
        if (x[k - 1] < x[k - 2]) return NAN;       /* must be sorted */

    if (ws) {
        W = ws;
    } else {
        w.mn  = malloc(((size_t)n + 1) * sizeof(int));
        w.mj  = malloc(((size_t)n + 1) * sizeof(int));
        w.gcm = malloc(((size_t)n + 1) * sizeof(int));
        w.lcm = malloc(((size_t)n + 1) * sizeof(int));
        if (!w.mn || !w.mj || !w.gcm || !w.lcm) {
            dip_ws_free(&w);
            return NAN;
        }
        W = &w;
    }
    /* 1-based views of the integer workspace (as in the Fortran original) */
    int *mn = W->mn - 1, *mj = W->mj - 1, *gcm = W->gcm - 1, *lcm = W->lcm - 1;
    const double *xp = x - 1;

    low = 1; high = n;
    /* N.B. dip_val is 2n * dip until the very end (Maechler's speedup). */
    dip_val = min_is_0 ? 0. : 1.;
    if (n < 2 || xp[n] == xp[1]) goto L_END;

    /* Indices over which combination is necessary for the convex MINORANT. */
    mn[1] = 1;
    for (j = 2; j <= n; ++j) {
        mn[j] = j - 1;
        for (;;) {
            mnj = mn[j];
            mnmnj = mn[mnj];
            if (mnj == 1 ||
                (xp[j] - xp[mnj]) * (mnj - mnmnj) <
                (xp[mnj] - xp[mnmnj]) * (j - mnj)) break;
            mn[j] = mnmnj;
        }
    }

    /* Indices over which combination is necessary for the concave MAJORANT. */
    mj[n] = n;
    for (k = n - 1; k >= 1; --k) {
        mj[k] = k + 1;
        for (;;) {
            mjk = mj[k];
            mjmjk = mj[mjk];
            if (mjk == n ||
                (xp[k] - xp[mjk]) * (mjk - mjmjk) <
                (xp[mjk] - xp[mjmjk]) * (k - mjk)) break;
            mj[k] = mjmjk;
        }
    }

LOOP_Start:
    /* Change points of the GCM from HIGH down to LOW. */
    gcm[1] = high;
    for (i = 1; gcm[i] > low; i++) gcm[i + 1] = mn[gcm[i]];
    ig = l_gcm = i;
    ix = ig - 1;

    /* Change points of the LCM from LOW up to HIGH. */
    lcm[1] = low;
    for (i = 1; lcm[i] < high; i++) lcm[i + 1] = mj[lcm[i]];
    ih = l_lcm = i;
    iv = 2;

    /* Largest distance > dip between the GCM and the LCM from LOW to HIGH. */
    {
        long double d = 0.;
        if (l_gcm != 2 || l_lcm != 2) {
            do {
                long double dx;
                int gcmix = gcm[ix], lcmiv = lcm[iv];
                if (gcmix > lcmiv) {
                    /* the next point of either the GCM or the LCM is from the LCM */
                    int gcmi1 = gcm[ix + 1];
                    dx = (lcmiv - gcmi1 + 1) -
                         ((long double)xp[lcmiv] - xp[gcmi1]) *
                             (gcmix - gcmi1) / (xp[gcmix] - xp[gcmi1]);
                    ++iv;
                    if (dx >= d) { d = dx; ig = ix + 1; ih = iv - 1; }
                } else {
                    /* ... or from the GCM.  Fix by Yong Lu (2003). */
                    int lcmiv1 = lcm[iv - 1];
                    dx = ((long double)xp[gcmix] - xp[lcmiv1]) * (lcmiv - lcmiv1) /
                             (xp[lcmiv] - xp[lcmiv1]) -
                         (gcmix - lcmiv1 - 1);
                    --ix;
                    if (dx >= d) { d = dx; ig = ix + 1; ih = iv; }
                }
                if (ix < 1) ix = 1;
                if (iv > l_lcm) iv = l_lcm;
            } while (gcm[ix] != lcm[iv]);
        } else {
            d = min_is_0 ? 0. : 1.;
        }
        if (d < dip_val) goto L_END;
    }

    /* Dips implied by the current LOW and HIGH. */
    {
        int j_best, j_l = -1, j_u = -1;

        dip_l = 0.;
        for (j = ig; j < l_gcm; ++j) {
            double max_t = 1.;
            int j_ = -1, jb = gcm[j + 1], je = gcm[j];
            if (je - jb > 1 && xp[je] != xp[jb]) {
                double C = (je - jb) / (xp[je] - xp[jb]);
                int jj;
                for (jj = jb; jj <= je; ++jj) {
                    double t = (jj - jb + 1) - (xp[jj] - xp[jb]) * C;
                    if (max_t < t) { max_t = t; j_ = jj; }
                }
            }
            if (dip_l < max_t) { dip_l = max_t; j_l = j_; }
        }

        dip_u = 0.;
        for (j = ih; j < l_lcm; ++j) {
            double max_t = 1.;
            int j_ = -1, jb = lcm[j], je = lcm[j + 1];
            if (je - jb > 1 && xp[je] != xp[jb]) {
                double C = (je - jb) / (xp[je] - xp[jb]);
                int jj;
                for (jj = jb; jj <= je; ++jj) {
                    double t = (xp[jj] - xp[jb]) * C - (jj - jb - 1);
                    if (max_t < t) { max_t = t; j_ = jj; }
                }
            }
            if (dip_u < max_t) { dip_u = max_t; j_u = j_; }
        }

        if (dip_u > dip_l) { dipnew = dip_u; j_best = j_u; }
        else               { dipnew = dip_l; j_best = j_l; }
        (void)j_best;
        if (dip_val < dipnew) dip_val = dipnew;
    }

    /* Necessary: otherwise infinite loop for unimodal samples (Maechler 1994). */
    if (low == gcm[ig] && high == lcm[ih]) goto L_END;
    low  = gcm[ig];
    high = lcm[ih];
    goto LOOP_Start;

L_END:
    if (lo) *lo = low;
    if (hi) *hi = high;
    if (!ws) dip_ws_free(&w);
    return (double)(dip_val / (2.0 * (double)n));
}

/* ================================================================== */
/* Per-pair accumulator                                                */
/* ================================================================== */

typedef struct {
    double *d;      /* distances, one slot per direction row */
    double *lr;     /* matching lr_ub (NAN when missing) */
    size_t n, cap;
} Vec;

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) { fprintf(stderr, "dip-pairs: out of memory\n"); exit(2); }
    return p;
}

static void vec_push(Vec *v, double d, double lr)
{
    if (v->n == v->cap) {
        size_t ncap = v->cap ? v->cap * 2 : 32;
        double *nd = realloc(v->d, ncap * sizeof(double));
        double *nl = realloc(v->lr, ncap * sizeof(double));
        if (!nd || !nl) { fprintf(stderr, "dip-pairs: out of memory\n"); exit(2); }
        v->d = nd; v->lr = nl; v->cap = ncap;
    }
    v->d[v->n] = d;
    v->lr[v->n] = lr;
    v->n++;
}

typedef struct {
    uint64_t h;
    char *config, *X, *Y;   /* X > Y */
    Vec xy, yx;
} Entry;

typedef struct { Entry *e; int64_t n, cap; } EntryArr;

static EntryArr g_out;

static void out_push(Entry *e)
{
    if (g_out.n == g_out.cap) {
        g_out.cap = g_out.cap ? g_out.cap * 2 : 1024;
        g_out.e = realloc(g_out.e, g_out.cap * sizeof(Entry));
        if (!g_out.e) { fprintf(stderr, "dip-pairs: out of memory\n"); exit(2); }
    }
    g_out.e[g_out.n++] = *e;
}

/* ================================================================== */
/* Hash table keyed by (config, X, Y), X > Y                           */
/* ================================================================== */

static Entry *tab;
static size_t tab_cap, tab_mask, tab_used;

static uint64_t fnv1a(const char *s)
{
    uint64_t h = 1469598103934665603ULL;
    while (*s) { h ^= (unsigned char)*s++; h *= 1099511628211ULL; }
    return h;
}

static uint64_t mix3(uint64_t a, uint64_t b, uint64_t c)
{
    uint64_t h = a;
    h ^= b + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h ^= c + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

/* Our own strdup(): avoids depending on the C library's, which ISO C does not
 * declare (and which would then only be visible via the feature-test macros
 * above). */
static char *xstrdup(const char *s)
{
    size_t n = strlen(s) + 1;
    char *p = malloc(n);
    if (!p) { fprintf(stderr, "dip-pairs: out of memory\n"); exit(2); }
    memcpy(p, s, n);
    return p;
}

static void tab_grow(void)
{
    size_t newcap = tab_cap ? tab_cap * 2 : (size_t)1 << 16;
    Entry *nt = calloc(newcap, sizeof(Entry));
    if (!nt) { fprintf(stderr, "dip-pairs: out of memory\n"); exit(2); }
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

static Entry *get_entry(const char *cfg, const char *X, const char *Y, int create)
{
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
    e->config = xstrdup(cfg);
    e->X = xstrdup(X);
    e->Y = xstrdup(Y);
    tab_used++;
    return e;
}

/* ================================================================== */
/* Per-pair statistics                                                 */
/* ================================================================== */

typedef struct {
    double dip, portion, max_unfiltered, max_selected;
    int64_t n_dip, n_kept, n_unfiltered, num_filtered, num_NA, num_lr_ub_zero;
    int used_filtered;
} Stats;

/* ascending by d; NaN sorts last */
typedef struct { double d, lr; } Side;

static int cmp_side(const void *a, const void *b)
{
    const Side *p = a, *q = b;
    int ax = isnan(p->d), ay = isnan(q->d);
    if (ax || ay) { if (ax && ay) return 0; return ax ? 1 : -1; }
    return (p->d < q->d) ? -1 : (p->d > q->d);
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    int ax = isnan(x), ay = isnan(y);
    if (ax || ay) { if (ax && ay) return 0; return ax ? 1 : -1; }
    return (x < y) ? -1 : (x > y);
}

static Side *side_sorted(const Vec *v)
{
    Side *s = xmalloc((v->n ? v->n : 1) * sizeof(Side));
    for (size_t i = 0; i < v->n; i++) { s[i].d = v->d[i]; s[i].lr = v->lr[i]; }
    qsort(s, v->n, sizeof(Side), cmp_side);
    return s;
}

static double bad2nan(const char *s)
{
    if (!s || !*s) return NAN;
    if (!strcmp(s, ".") || !strcmp(s, "nan") || !strcmp(s, "NaN")
        || !strcmp(s, "NA") || !strcmp(s, "-")) return NAN;
    char *end;
    double v = strtod(s, &end);
    if (end == s || *end) return NAN;
    return v;
}

static FILE *g_dump;   /* --dump-selected sink (main thread only!) */

/* Reconcile the two directions, apply the lr_ub filter and the portion switch.
 * On return *sel points at the selected distances (ascending, no NaN) and
 * *nsel is their count; *sel must be freed by the caller.  All other fields of
 * st except dip are filled in. */
static double *select_distances(const Entry *e, double chi_sq, double min_portion,
                                Stats *st, size_t *nsel_out)
{
    double *sel = NULL;
    *nsel_out = 0;
    st->dip = NAN; st->portion = NAN;
    st->max_unfiltered = NAN; st->max_selected = NAN;

    size_t n1 = e->xy.n, n2 = e->yx.n;
    size_t nmax = n1 > n2 ? n1 : n2;
    if (nmax == 0) return NULL;

    Side *s1 = side_sorted(&e->xy);
    Side *s2 = side_sorted(&e->yx);

    double *full = xmalloc(nmax * sizeof(double));   /* unfiltered reconciled d */
    double *kept = xmalloc(nmax * sizeof(double));   /* filtered reconciled d   */
    size_t nf = 0, nk = 0;

    for (size_t i = 0; i < nmax; i++) {
        int has1 = i < n1 && !isnan(s1[i].d);
        int has2 = i < n2 && !isnan(s2[i].d);
        const Side *c;
        if (has1 && has2) c = (s1[i].d <= s2[i].d) ? &s1[i] : &s2[i];
        else if (has1)    c = &s1[i];
        else if (has2)    c = &s2[i];
        else { st->num_NA++; continue; }             /* reconciled row is NA */
        if (isnan(c->d)) { st->num_NA++; continue; }

        if (c->lr == 0.0) st->num_lr_ub_zero++;
        full[nf++] = c->d;
        if (isnan(st->max_unfiltered) || c->d > st->max_unfiltered)
            st->max_unfiltered = c->d;

        /* keep strictly lr_ub > chi_sq; missing/NA lr_ub is dropped */
        if (!isnan(c->lr) && c->lr > chi_sq) kept[nk++] = c->d;
    }
    free(s1); free(s2);

    st->n_unfiltered = (int64_t)nf;
    st->n_kept = (int64_t)nk;
    st->num_filtered = (int64_t)(nf - nk);
    st->portion = nf ? (double)nk / (double)nf : NAN;

    /* the iddd switch: filtered only if the kept portion exceeds min_portion */
    int use_filtered = (nf > 0) && (st->portion > min_portion);
    st->used_filtered = use_filtered;
    *nsel_out = use_filtered ? nk : nf;
    st->n_dip = (int64_t)*nsel_out;

    /* the selected set is ascending in either case (both inputs are sorted) */
    sel = xmalloc((*nsel_out ? *nsel_out : 1) * sizeof(double));
    memcpy(sel, use_filtered ? kept : full, *nsel_out * sizeof(double));
    for (size_t i = 0; i < *nsel_out; i++)
        if (isnan(st->max_selected) || sel[i] > st->max_selected)
            st->max_selected = sel[i];

    free(full); free(kept);
    return sel;
}

/* Reconcile + filter + dip for one pair.  Safe to call from several threads:
 * it touches no shared state (the --dump-selected sink is written separately). */
static Stats compute_pair(const Entry *e, double chi_sq, double min_portion,
                          int64_t min_n)
{
    Stats st;
    memset(&st, 0, sizeof(st));
    size_t nsel = 0;
    double *sel = select_distances(e, chi_sq, min_portion, &st, &nsel);
    if (!sel) return st;

    if ((int64_t)nsel >= min_n && nsel >= 2) {
        /* dip() indexes with int; ignore absurdly large single pairs */
        if (nsel > (size_t)INT32_MAX) {
            fprintf(stderr, "dip-pairs: warning: pair %s/%s has %zu selected rows; "
                            "capping at %d for the dip statistic\n",
                    e->X, e->Y, nsel, INT32_MAX);
            nsel = (size_t)INT32_MAX;
        }
        int n = (int)nsel;
        double *xs = xmalloc(((size_t)n + 1) * sizeof(double));
        int m = 0;
        for (size_t i = 0; i < nsel; i++) if (!isnan(sel[i])) xs[m++] = sel[i];
        if (m >= 2 && (int64_t)m >= min_n) {
            /* sel is already ascending; guard against a caller breaking that */
            for (int i = 1; i < m; i++)
                if (xs[i] < xs[i - 1]) {
                    qsort(xs, (size_t)m, sizeof(double), cmp_double);
                    break;
                }
            st.dip = dip(xs, m, NULL, NULL, NULL, 0);
            st.n_dip = m;
        } else {
            st.n_dip = m;
        }
        free(xs);
    }

    free(sel);
    return st;
}

/* Write one --dump-selected row for one pair.  Must run on the main thread:
 * the parallel stats pass deliberately does not touch g_dump. */
static void dump_selected(const Entry *e, double chi_sq, double min_portion)
{
    Stats st;
    size_t nsel = 0;
    memset(&st, 0, sizeof(st));
    double *sel = select_distances(e, chi_sq, min_portion, &st, &nsel);
    fprintf(g_dump, "%s\t%s\t%s\t", e->config, e->X, e->Y);
    for (size_t i = 0; i < nsel; i++)
        fprintf(g_dump, "%s%.17g", i ? "," : "", sel[i]);
    fputs("\n", g_dump);
    free(sel);
}

/* ================================================================== */
/* Output                                                              */
/* ================================================================== */

#define OUT_HDR \
    "config\tgenome_a\tgenome_b\tdip\tn_dip\tportion\tn_kept\tn_unfiltered" \
    "\tused\tnum_filtered\tnum_NA\tmax_unfiltered_distance" \
    "\tmax_dip_distance\tnum_lr_ub_zero\n"

static int g_prec = 10;   /* significant digits for floating point output */

static void emit(const Entry *e, const Stats *st)
{
    char fmt[8];
    printf("%s\t%s\t%s\t", e->config, e->X, e->Y);
    if (isnan(st->dip)) printf("nan"); else printf("%.*g", g_prec, st->dip);
    printf("\t%lld\t", (long long)st->n_dip);
    if (isnan(st->portion)) printf("nan"); else printf("%.*g", g_prec, st->portion);
    printf("\t%lld\t%lld\t%s\t%lld\t%lld\t",
           (long long)st->n_kept, (long long)st->n_unfiltered,
           st->used_filtered ? "filtered" : "unfiltered",
           (long long)st->num_filtered, (long long)st->num_NA);
    snprintf(fmt, sizeof(fmt), "%%.%dg", g_prec);
    if (isnan(st->max_unfiltered)) printf("nan"); else printf(fmt, st->max_unfiltered);
    printf("\t");
    if (isnan(st->max_selected)) printf("nan"); else printf(fmt, st->max_selected);
    printf("\t%lld\n", (long long)st->num_lr_ub_zero);
}

/* ================================================================== */
/* Input parsing                                                       */
/* ================================================================== */

static int split(char *line, char **f, int maxf)
{
    int n = 0;
    char *p = line;
    while (n < maxf) {
        char *t = strchr(p, '\t');
        if (!t) { f[n++] = p; break; }
        *t = 0; f[n++] = p; p = t + 1;
    }
    return n;
}

static int cmp_entry(const void *a, const void *b)
{
    const Entry *e1 = a, *e2 = b;
    int r;
    if ((r = strcmp(e1->config, e2->config))) return r;
    if ((r = strcmp(e1->X, e2->X))) return r;
    return strcmp(e1->Y, e2->Y);
}

/* stream mode: rows of one pair should be contiguous.  The last RECENT flushed
 * raw keys are remembered so that a late row of an already-emitted pair is
 * reported instead of silently duplicating it. */
#define RECENT 1024
typedef struct { uint64_t h; char *cfg, *a, *b; } RawKey;
static RawKey recent[RECENT];
static size_t recent_next;
static int64_t n_dupe_rows;

static int raw_seen(uint64_t h, const char *cfg, const char *a, const char *b)
{
    for (size_t i = 0; i < RECENT; i++)
        if (recent[i].h == h && recent[i].cfg &&
            !strcmp(recent[i].cfg, cfg) && !strcmp(recent[i].a, a)
            && !strcmp(recent[i].b, b)) return 1;
    return 0;
}

static void raw_remember(uint64_t h, const char *cfg, const char *a, const char *b)
{
    RawKey *r = &recent[recent_next];
    recent_next = (recent_next + 1) % RECENT;
    free(r->cfg); free(r->a); free(r->b);
    r->h = h; r->cfg = xstrdup(cfg); r->a = xstrdup(a); r->b = xstrdup(b);
}

static void usage(const char *argv0)
{
    fprintf(stderr,
        "usage: %s <gdiff-sample.tsv> [--chi-sq 3.841] [--min-portion 0.66]\n"
        "       [--min-n 4] [--jobs N] [--stream|--no-stream] [--dump-selected F]\n"
        "       [--precision N]\n"
        "\n"
        "Hartigan's dip statistic for every symmetric pair of a gdiff sample file,\n"
        "after the same reconcile/filter/selection strategy as estimate-sample-ani.\n"
        "See the header comment of this file for the exact semantics.\n", argv0);
}

#ifdef DIP_PAIRS_LIB
int dip_pairs_main(int argc, char **argv)
#else
int main(int argc, char **argv)
#endif
{
    if (argc < 2) { usage(argv[0]); return 1; }
    const char *path = argv[1];
    const char *dump_path = NULL;
    double chi_sq = 3.841, min_portion = 0.66;
    int64_t min_n = 2;
    int jobs = 0, stream = -1;   /* -1: auto by input size */

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--chi-sq") && i + 1 < argc) chi_sq = atof(argv[++i]);
        else if (!strcmp(argv[i], "--min-portion") && i + 1 < argc) min_portion = atof(argv[++i]);
        else if (!strcmp(argv[i], "--min-n") && i + 1 < argc) min_n = atoll(argv[++i]);
        else if (!strcmp(argv[i], "--jobs") && i + 1 < argc) jobs = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--stream")) stream = 1;
        else if (!strcmp(argv[i], "--no-stream")) stream = 0;
        else if (!strcmp(argv[i], "--dump-selected") && i + 1 < argc) dump_path = argv[++i];
        else if (!strcmp(argv[i], "--precision") && i + 1 < argc) g_prec = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(argv[0]); return 0; }
        else { fprintf(stderr, "dip-pairs: unknown arg: %s\n", argv[i]); usage(argv[0]); return 1; }
    }
    if (jobs < 0) jobs = 1;
#ifdef _OPENMP
    if (jobs == 0) jobs = omp_get_max_threads();
    if (jobs > 1) omp_set_num_threads(jobs);
#else
    if (jobs == 0) jobs = 1;
    jobs = 1;
#endif

    struct stat sb;
    if (stream == -1)
        stream = (stat(path, &sb) == 0 && sb.st_size >= STREAM_THRESHOLD);

    FILE *fh = fopen(path, "r");
    if (!fh) { fprintf(stderr, "dip-pairs: cannot open %s\n", path); return 1; }
    if (dump_path) {
        g_dump = fopen(dump_path, "w");
        if (!g_dump) { fprintf(stderr, "dip-pairs: cannot write %s\n", dump_path); return 1; }
    }

    tab_grow();

    if (stream) { printf(OUT_HDR); fflush(stdout); }

    char *line = NULL; size_t linelen = 0; ssize_t nrl;
    size_t lineno = 0;
    int64_t n_rows = 0, n_skipped = 0, n_selfpairs = 0;
    Entry *cur = NULL;                 /* stream mode: pair being accumulated */
    int64_t n_stream_pairs = 0;

    while ((nrl = getline(&line, &linelen, fh)) != -1) {
        lineno++;
        {
            /* skip a header row if present (first line starting with "config\t") */
            const char *h = "config\t";
            if (lineno == 1 && nrl >= 7 && !strncmp(line, h, 7)) continue;
        }
        while (nrl > 0 && (line[nrl - 1] == '\n' || line[nrl - 1] == '\r'))
            line[--nrl] = '\0';
        if (nrl == 0) continue;
        n_rows++;
        char *f[11];
        int nf = split(line, f, 11);
        if (nf < 9) { n_skipped++; continue; }
        char *cfg = f[0], *ga = f[1], *gb = f[2];
        int c = strcmp(ga, gb);
        if (c == 0) { n_selfpairs++; continue; }
        const char *X = c > 0 ? ga : gb;
        const char *Y = c > 0 ? gb : ga;
        double d = bad2nan(f[8]);
        double lr = bad2nan(nf > 10 ? f[10] : NULL);

        if (stream) {
            uint64_t rh = mix3(fnv1a(cfg), fnv1a(ga), fnv1a(gb));
            if (cur && !strcmp(cur->config, cfg) && !strcmp(cur->X, X)
                && !strcmp(cur->Y, Y)) {
                /* the current pair continues */
            } else if (raw_seen(rh, cfg, ga, gb)) {
                n_dupe_rows++;         /* late row of an already emitted pair */
                continue;
            } else {
                if (cur) {
                    Stats st = compute_pair(cur, chi_sq, min_portion, min_n);
                    emit(cur, &st);
                    if (g_dump) dump_selected(cur, chi_sq, min_portion);
                    n_stream_pairs++;
                    raw_remember(cur->h, cur->config, cur->X, cur->Y);
                    /* also remember (Y,X) so a trailing reversed row is caught */
                    raw_remember(mix3(fnv1a(cur->config), fnv1a(cur->Y), fnv1a(cur->X)),
                                 cur->config, cur->Y, cur->X);
                    cur->xy.n = cur->yx.n = 0;      /* keep the buffers */
                } else {
                    cur = xmalloc(sizeof(Entry));
                    memset(cur, 0, sizeof(Entry));
                }
                free(cur->config); free(cur->X); free(cur->Y);
                cur->config = xstrdup(cfg);
                cur->X = xstrdup(X);
                cur->Y = xstrdup(Y);
                cur->h = mix3(fnv1a(cfg), fnv1a(X), fnv1a(Y));
            }
            vec_push(c > 0 ? &cur->xy : &cur->yx, d, lr);
        } else {
            Entry *e = get_entry(cfg, X, Y, 1);
            vec_push(c > 0 ? &e->xy : &e->yx, d, lr);
        }
    }
    free(line);
    fclose(fh);

    if (stream) {
        if (cur) {
            Stats st = compute_pair(cur, chi_sq, min_portion, min_n);
            emit(cur, &st);
            if (g_dump) dump_selected(cur, chi_sq, min_portion);
            n_stream_pairs++;
        }
        fflush(stdout);
        if (n_dupe_rows)
            fprintf(stderr, "dip-pairs: warning: %lld row(s) arrived after their pair "
                            "was emitted; input is not grouped by pair -- use "
                            "--no-stream\n", (long long)n_dupe_rows);
    } else {
        for (size_t i = 0; i < tab_cap; i++)
            if (tab[i].config) out_push(&tab[i]);
        if (g_out.n > 1) qsort(g_out.e, (size_t)g_out.n, sizeof(Entry), cmp_entry);

        int64_t n = g_out.n;
        Stats *sts = xmalloc((size_t)(n ? n : 1) * sizeof(Stats));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) if(jobs > 1 && n > 1)
#endif
        for (int64_t i = 0; i < n; i++)
            sts[i] = compute_pair(&g_out.e[i], chi_sq, min_portion, min_n);

        printf(OUT_HDR);
        for (int64_t i = 0; i < n; i++) {
            emit(&g_out.e[i], &sts[i]);
            if (g_dump) dump_selected(&g_out.e[i], chi_sq, min_portion);
        }
        fflush(stdout);
        free(sts);
    }

    if (g_dump) fclose(g_dump);
    fprintf(stderr,
            "dip-pairs: %lld pairs from %lld rows (%lld skipped, %lld self-pairs, "
            "%lld late) -- %s mode, %d thread(s)\n",
            (long long)(stream ? n_stream_pairs : g_out.n), (long long)n_rows,
            (long long)n_skipped, (long long)n_selfpairs, (long long)n_dupe_rows,
            stream ? "stream" : "vector", jobs);
    return 0;
}
