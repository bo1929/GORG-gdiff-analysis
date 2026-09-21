# Resource benchmarking

All-by-all pairwise distance benchmark over a random sample of genomes, one
config per method, measuring wall time, CPU time and peak RSS.

Edit the defaults at the top of `run_benchmarks.sh` (threads, sample size,
per-method args), then:

```bash
cd resource-benchmarking
./run_benchmarks.sh
```

Everything lands in `output/`: the genome sample, sketches, distance files,
`/usr/bin/time` reports, and `output/resources.tsv` (one row per phase).

Run a subset with `./run_benchmarks.sh dashing2 skani`. Set `FORCE=1` at the
top of the script to redraw the sample and redo methods already in
`resources.tsv`.

Default methods: `dashing2 skani mash gdiff fastani`. `fastani` must be on
`PATH` (`fastani` or `fastANI`). `anib` is implemented; add it to `METHODS`
to run it.

## `output/resources.tsv`

| column | meaning |
|---|---|
| `method` | `dashing2`, `skani`, `mash`, `gdiff`, `fastani`, `anib` |
| `config` | preset used (e.g. `v4`, `slow`, `abcs`) |
| `phase` | `sketch` or `dist` |
| `threads` | value of `THREADS` |
| `wall_sec` | elapsed wall-clock seconds |
| `user_sec` | user CPU seconds |
| `sys_sec` | system CPU seconds |
| `maxrss_mib` | peak RSS of the process tree, MiB |

| method | sketch | dist |
|---|---|---|
| dashing2 | `dashing2 sketch --cache-sketches` | `dashing2 cmp` over the cached sketches |
| skani | `skani sketch --separate-sketches` | `skani triangle -l <sketches>` |
| gdiff | `gdiff sketch` one bundle | `gdiff dist <bundle>` |
| mash | `mash sketch -l` one `.msh` | `mash dist all.msh all.msh` |
| fastani | *(none)* | `fastani --ql --rl` (sketches internally) |
| anib | *(none)* | `anib <list>` |

Shipped presets: dashing2 `v4` (`--symmetric-containment -k 23 -S 2048`),
skani `slow` (`--slow`), gdiff `abcs` (`-k 23 -w 23 --frac 0.5 -l 500
--sample-size 1000` / `--hdist-th 3`), mash `default` (`-k 21 -s 1000`).
Change them in the defaults block.

Wall/CPU/RSS come from `/usr/bin/time` (`-l` on macOS, `-v` on Linux).
Distance files are written so the run is verifiable; only resource usage is
benchmarked. `skani triangle` writes a lower-triangular PHYLIP matrix (plus
a sibling `.af` file). `mash dist msh msh` writes the full `N x N` table.
`fastani --ql --rl` writes `N x N` including self comparisons.
