/* dip_only.c -- compute Hartigan's dip for vectors dumped by dip-pairs.
 *
 * Input: the --dump-selected format of scripts/dip-pairs, one line per pair:
 *     config <TAB> X <TAB> Y <TAB> d,d,d,...
 * Output: config <TAB> X <TAB> Y <TAB> n <TAB> dip
 *
 * This exists to cross-check the dip implementation independently of the
 * reconcile/filter stage (see scripts/validate_dip_vectors.sh) and is a handy
 * way to run the statistic over pre-extracted samples.
 *
 * Build:  cc -O3 -o scripts/dip-only scripts/dip_only.c -lm
 * Usage:  dip-only <selected.tsv> [--precision N]
 */
#define DIP_PAIRS_LIB 1
#include "dip_pairs.c"

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s <dump-selected.tsv> [--precision N]\n", argv[0]); return 1; }
    int prec = 10;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--precision") && i + 1 < argc) prec = atoi(argv[++i]);
        else { fprintf(stderr, "dip-only: unknown arg %s\n", argv[i]); return 1; }
    }
    FILE *fh = fopen(argv[1], "r");
    if (!fh) { fprintf(stderr, "dip-only: cannot open %s\n", argv[1]); return 1; }

    char *line = NULL; size_t cap = 0; ssize_t len;
    printf("config\tgenome_a\tgenome_b\tn\tdip\n");
    while ((len = getline(&line, &cap, fh)) != -1) {
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;
        char *f[4];
        if (split(line, f, 4) < 4) continue;
        size_t n = 0;
        double *v = NULL;
        if (*f[3]) {
            size_t a_cap = 64;
            v = xmalloc(a_cap * sizeof(double));
            char *p = f[3];
            while (*p) {
                if (n == a_cap) { a_cap *= 2; v = realloc(v, a_cap * sizeof(double)); if (!v) return 2; }
                char *end = p;
                double val = strtod(p, &end);
                if (end == p) {
                    /* not a number: skip to the next comma (never spin) */
                    fprintf(stderr, "dip-only: skipping non-numeric field at "
                                    "offset %ld: %.16s\n", (long)(p - f[3]), p);
                    while (*p && *p != ',') p++;
                    if (*p == ',') p++;
                    continue;
                }
                p = end;
                v[n++] = val;
                if (*p == ',') p++;
            }
        }
        double d = (n >= 2) ? dip(v, (int)n, NULL, NULL, NULL, 0) : NAN;
        printf("%s\t%s\t%s\t%zu\t", f[0], f[1], f[2], n);
        if (isnan(d)) printf("nan"); else printf("%.*g", prec, d);
        printf("\n");
        free(v);
    }
    free(line);
    fclose(fh);
    return 0;
}
