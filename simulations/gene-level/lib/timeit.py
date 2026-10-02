#!/usr/bin/env python3
"""Run a command under /usr/bin/time and append one resource row to a TSV.

    python3 timeit.py --tsv resources.tsv --pair VARIANT --label wfmash --threads 1 -- wfmash ...

Appends::

    pair  label  threads  wall_s  user_s  sys_s  max_rss_mb  exit

``--pair`` identifies which simulated pair the row belongs to; without it several pairs'
rows merge into one indistinguishable table, so the batch drivers always pass it.

macOS ``time -l`` reports maximum resident set size in bytes and GNU ``time -v`` in kbytes;
both are normalised to MiB here. The command's own stdout/stderr are passed through, so the
usual redirections still apply. The wrapped command's exit status becomes this script's.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile

DARWIN_RE = re.compile(r"^\s*([\d.]+)\s+real\s+([\d.]+)\s+user\s+([\d.]+)\s+sys")
RSS_RE = re.compile(r"^\s*(\d+)\s+maximum resident set size", re.M)
GNU_RE = re.compile(
    r"User time \(seconds\):\s*([\d.]+).*?"
    r"System time \(seconds\):\s*([\d.]+).*?"
    r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+).*?"
    r"Maximum resident set size \(kbytes\):\s*(\d+)",
    re.S,
)
FIELDS = ["pair", "label", "threads", "wall_s", "user_s", "sys_s", "max_rss_mb", "exit"]


def hmstime(tok: str) -> float:
    parts = [float(p) for p in tok.split(":")]
    s = 0.0
    for p in parts:
        s = s * 60 + p
    return s


def parse(text: str):
    """-> (wall_s, user_s, sys_s, max_rss_mb) or None."""
    m = DARWIN_RE.match(text)
    if m:
        wall, user, sysd = (float(m.group(i)) for i in (1, 2, 3))
        r = RSS_RE.search(text)
        return wall, user, sysd, (int(r.group(1)) / 2**20 if r else float("nan"))
    m = GNU_RE.search(text)
    if m:
        user, sysd = float(m.group(1)), float(m.group(2))
        wall = hmstime(m.group(3))
        return wall, user, sysd, int(m.group(4)) / 1024
    return None


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tsv", required=True, help="TSV to append to (header written if absent)")
    p.add_argument("--label", required=True, help="row label, e.g. wfmash")
    p.add_argument("--pair", default="-", help="pair/variant this row belongs to [-]")
    p.add_argument("--threads", default="1", help="threads the command was told to use [1]")
    p.add_argument("cmd", nargs=argparse.REMAINDER)
    args = p.parse_args()

    cmd = args.cmd[1:] if args.cmd and args.cmd[0] == "--" else args.cmd
    if not cmd:
        raise SystemExit("no command given; use: timeit.py ... -- <cmd> [args]")

    flag = "-l" if sys.platform == "darwin" else "-v"
    with tempfile.NamedTemporaryFile("r", suffix=".time", delete=False) as tf:
        report = tf.name
    try:
        rc = subprocess.call(["/usr/bin/time", flag, "-o", report] + cmd)
        text = open(report).read()
    finally:
        os.unlink(report)

    row = parse(text)
    if row is None:
        print(f"timeit: could not parse the time report for {args.label}:\n{text}",
              file=sys.stderr)
        row = (float("nan"),) * 4
    wall, user, sysd, rss = row

    new = not os.path.exists(args.tsv) or os.path.getsize(args.tsv) == 0
    with open(args.tsv, "a") as fh:
        if new:
            fh.write("\t".join(FIELDS) + "\n")
        fh.write(f"{args.pair}\t{args.label}\t{args.threads}\t{wall:.3f}\t{user:.3f}\t"
                 f"{sysd:.3f}\t{rss:.1f}\t{rc}\n")
    print(f"    {args.label}: {wall:.2f}s wall, {rss:.1f} MiB peak, exit {rc}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
