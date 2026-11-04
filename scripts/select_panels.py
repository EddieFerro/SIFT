#!/usr/bin/env python3
"""SIFT panel selection: rank the pangenome's haplotypes for one assembly from the pseudo-MEMs of its sequences and
write donor panels (e.g. for ImpuT2T). Needs only the output of `sift query`; no ground truth.

    python3 select_panels.py hits/chr*.L1000.hits.gz --contigs contigs/ --exclude HG002 \
        --always CHM13 GRCh38 --n 10 25 50 --mode both --out HG002.select

Input: one or more files with the FASTA headers of `sift query` (gzipped or not), e.g.
    sift query contigs/chr1.fa -d panel/chr1 -i index/chr1 -L 1000 -o /dev/stdout | grep '^>' | gzip > hits/chr1.L1000.hits.gz
Header: >NAME:START-END tags=sample#hap,...   (other lines are ignored)
  Groups: each file is one group, ranked separately. A file whose name contains a chromosome (chr1 ... chrY, chrM) is
          in the group of that chromosome; otherwise the group is the file name without its extensions (.gz, .hits,
          .L<n>, .txt). Files of the same group are pooled. With a single whole-genome hits file, the one group's
          ranking is also the genome-wide ranking.
  Names:  a record NAME of the form  query_hap|contig|+  or  query_hap|contig|-  (as the pipeline's
          extract_contigs.py writes: both strands of a contig as separate records) is scored as one contig, the -
          strand's MEMs mapped to forward coordinates. Any other NAME is one forward-strand sequence.
  --contigs gives exact sequence lengths: FASTA files, or folders holding <group>.fa[.gz] / <group>.fasta[.gz].
          Without it, a sequence's length is estimated by its last MEM end (affects only the contig-end windows).

Ranking, per group:
  cov           bp of the group's sequences covered by pseudo-MEMs tagged with the haplotype (both strands)
  ends_idf      the same in the --window bp at each sequence end, each MEM's bp weighted log(H/|tags|) (H haplotypes),
                ranked per sequence end and combined over the group's sequence ends by reciprocal-rank fusion
  ranking       reciprocal-rank fusion (k = --rrf-k) of cov and ends_idf
Genome-wide ranking: reciprocal-rank fusion of the per-group rankings (the groups in --groups, default all).

Retained and excluded haplotypes (--always/--always-file, --exclude/--exclude-file): a name is a sample (all of its
haplotypes, e.g. HG002) or one haplotype (HG002.1, or PanSN HG002#1). Excluded haplotypes (the query, its relatives)
are dropped. Retained ones (e.g. the references CHM13 GRCh38) are not ranked and head every panel, in the order given.

Panels (each = the retained haplotypes + the top N - (number retained) ranked haplotypes):
  --mode chr     a panel per group, from that group's ranking                (SIFT-CHR)
  --mode genome  one panel for every group, from the genome-wide ranking     (SIFT-Genome)
  --mode both    both

Outputs (<out> = --out):
  <out>.rankings.tsv            variant chrom hap rank cov_bp   (variants: fuse200k = the ranking, chr_cov = cov;
                                the chrom column holds the group)
  <out>.genome_ranking.tsv      rank hap score
  <out>.sift_chr_N<N>.tsv       group n haplotypes(comma-separated)
  <out>.sift_genome_N<N>.txt    one haplotype per line
"""
import argparse, gzip, heapq, math, os, re, sys
from collections import defaultdict

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("hits", nargs="+", help="sift query headers, one or more files (see Groups)")
ap.add_argument("--contigs", nargs="*", default=[], help="FASTA files or folders of <group>.fa: exact sequence lengths")
ap.add_argument("--exclude", nargs="*", default=[], help="samples or haplotypes never ranked or retained")
ap.add_argument("--exclude-file", help="file of samples or haplotypes to exclude, one per line")
ap.add_argument("--always", nargs="*", default=[], help="samples or haplotypes in every panel, never ranked")
ap.add_argument("--always-file", help="file of samples or haplotypes to retain, one per line")
ap.add_argument("--window", type=int, default=200000, help="sequence-end window (bp)")
ap.add_argument("--rrf-k", type=float, default=10.0)
ap.add_argument("--groups", "--chroms", dest="groups", nargs="*",
                help="groups in the genome-wide ranking (default: all with hits)")
ap.add_argument("--n", type=int, nargs="*", default=[], help="panel sizes, the retained haplotypes included")
ap.add_argument("--mode", choices=["chr", "genome", "both"], default="both")
ap.add_argument("--out", required=True)
a = ap.parse_args()

def name_of(tag):                                   # PanSN sample#hap[#contig] -> sample (hap 0) or sample.hap
    if "#" not in tag: return tag
    s, h = tag.split("#")[:2]
    return s if h == "0" else f"{s}.{h}"
def sample(h): return h.rsplit(".", 1)[0] if "." in h else h
def read_list(path):
    if not path: return []
    out = []
    for line in open(path):
        f = line.split()
        if f and not f[0].startswith("#"): out.append(f[0])     # lines starting with # are comments
    return out
def names(xs): return [name_of(x) for x in xs]
always = list(dict.fromkeys(names(a.always + read_list(a.always_file))))
excl = set(names(a.exclude + read_list(a.exclude_file)))
both = set(always) & excl
if both: sys.exit(f"error: both retained and excluded: {' '.join(sorted(both))}")
alw = set(always)
def is_always(h): return h in alw or sample(h) in alw
def is_excluded(h): return h in excl or sample(h) in excl

CHROM = re.compile(r"(?:^|[._/-])(chr(?:[0-9]+|X|Y|M))(?=[._]|$)")
SUFFIX = re.compile(r"(\.(gz|bgz|hits|txt|tsv|fa|fasta|out)|\.L\d+)$")
def group_of(path):
    base = os.path.basename(path)
    m = CHROM.search(base)
    if m: return m.group(1)
    while True:
        b = SUFFIX.sub("", base)
        if b == base or not b: return base
        base = b
STRANDED = re.compile(r"^(.+)\|([^|]+)\|([+-])$")
def parse_name(name):                              # -> (query hap, sequence, strand)
    m = STRANDED.match(name)
    return m.groups() if m else ("", name, "+")
def opener(path): return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)

# ---- sequence lengths and hits ------------------------------------------------------------------------------
clen = {}
def read_lengths(path):
    key = None; n = 0
    for line in opener(path):
        if line[0] == ">":
            if key: clen[key] = n
            q, c, st = parse_name(line[1:].split()[0])
            key = (q, c) if st == "+" else None; n = 0
        else:
            n += len(line.strip())
    if key: clen[key] = n
FA = (".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz")
contig_dirs = [p for p in a.contigs if os.path.isdir(p)]
for p in a.contigs:
    if not os.path.isdir(p):
        if not os.path.exists(p): sys.exit(f"error: no such file or folder: {p}")
        read_lengths(p)
frags = defaultdict(list)
panel, seen = set(), set()
done_len = set()
for path in a.hits:
    grp = group_of(path)
    if contig_dirs and grp not in done_len:
        base = os.path.basename(path); i = base.find(grp)
        cands = ([base[:i] + grp + ext for ext in FA] if i > 0 else []) + [grp + ext for ext in FA]
        f = next((os.path.join(d, nm) for d in contig_dirs for nm in cands if os.path.exists(os.path.join(d, nm))), None)
        if f: read_lengths(f)
        done_len.add(grp)
    for line in opener(path):
        if line[0] != ">": continue
        f = line[1:].split()
        if not f or ":" not in f[0]: continue
        name, rng = f[0].rsplit(":", 1)
        se = rng.split("-")
        if len(se) != 2 or not se[0].isdigit() or not se[1].isdigit(): continue
        s, e = int(se[0]), int(se[1])
        q, c, st = parse_name(name)
        tags = next((x[5:] for x in f[1:] if x.startswith("tags=")), "")
        if not tags:
            clen.setdefault((q, c), e)
            continue
        hs = {name_of(t) for t in tags.split(",")}
        seen.update(hs)
        hs = tuple(sorted(h for h in hs if not is_always(h) and not is_excluded(h)))
        if not hs: continue
        frags[(grp, q, c)].append((s, e, st, hs))
        panel.update(hs)
H = len(panel)
if not H: sys.exit("error: no haplotypes to rank (no tagged hits left after --exclude/--always)")
est = 0
for (grp, q, c), fl in frags.items():
    if (q, c) not in clen:
        clen[(q, c)] = max(e for _, e, _, _ in fl); est += 1
if est: print(f"[select] warning: {est} sequence lengths estimated from hits (pass --contigs for exact)", file=sys.stderr)
mems = {}
for (grp, q, c), fl in frags.items():
    n = clen[(q, c)]
    mems[(grp, q, c)] = [(s - 1, e, st, hs) if st == "+" else (n - e, n - s + 1, st, hs) for s, e, st, hs in fl]

# retained haplotypes, in the order given: a sample expands to its haplotypes seen in the hits
REFS = []
for x in always:
    hit = sorted(h for h in seen if h == x or sample(h) == x)
    if not hit: print(f"[select] warning: retained {x} has no hits; written to the panels as given", file=sys.stderr)
    for h in hit or [x]:
        if h not in REFS: REFS.append(h)

# ---- the two scores fuse200k needs (cov, idf of rescore_hits.score_region) ------------------------
def union_len(xs):
    xs.sort(); tot = 0; cs = ce = -1
    for s, e in xs:
        if s > ce: tot += max(0, ce - cs); cs, ce = s, e
        else: ce = max(ce, e)
    return tot + max(0, ce - cs)
def weighted_union(xs):
    ev = sorted([(s, 0, w, i) for i, (s, e, w) in enumerate(xs)] + [(e, 1, w, i) for i, (s, e, w) in enumerate(xs)])
    heap, dead, tot, last = [], set(), 0.0, None
    for pos, kind, w, i in ev:
        while heap and heap[0][1] in dead: heapq.heappop(heap)
        if heap and last is not None: tot += -heap[0][0] * (pos - last)
        last = pos
        if kind == 0: heapq.heappush(heap, (-w, i))
        else: dead.add(i)
    return tot
def score_region(ms, lo, hi, want):
    per = defaultdict(lambda: {"all": [], "w": []})
    for s, e, st, hs in ms:
        s2, e2 = max(s, lo), min(e, hi)
        if e2 <= s2: continue
        w = math.log(H / len(hs)) if len(hs) < H else 0.0
        for h in hs:
            d = per[h]; d["all"].append((s2, e2)); d["w"].append((s2, e2, w))
    if want == "cov": return {h: union_len(list(d["all"])) for h, d in per.items()}
    return {h: weighted_union(d["w"]) for h, d in per.items()}
def ranks(score):
    items = sorted(((v, h) for h, v in score.items() if v > 0), reverse=True)
    out, prev, r = {}, None, 0
    for i, (v, h) in enumerate(items):
        if v != prev: r = i + 1; prev = v
        out[h] = r
    return out

# ---- per-group rankings ---------------------------------------------------------------------------
W, K = a.window, a.rrf_k
chroms = sorted({k[0] for k in mems})
chr_cov, fused_rank = {}, {}
for chrom in chroms:
    tot = defaultdict(float); ends = defaultdict(float)
    for (ch, q, c), ms in mems.items():
        if ch != chrom: continue
        n = clen[(q, c)]
        for h, v in score_region(ms, 0, n, "cov").items(): tot[h] += v
        for lo, hi in ((0, min(W, n)), (max(0, n - W), n)):
            for h, rk in ranks(score_region(ms, lo, hi, "idf")).items(): ends[h] += 1.0 / (K + rk)
    chr_cov[chrom] = dict(tot)
    base, er = ranks(tot), ranks(ends)
    fused = defaultdict(float)
    for rk in (base, er):
        for h, r in rk.items(): fused[h] += 1.0 / (K + r)
    fused_rank[chrom] = (ranks(fused), base)
    print(f"[select] {chrom}: {len(fused_rank[chrom][0])} haplotypes ranked", file=sys.stderr)

with open(a.out + ".rankings.tsv", "w") as f:
    f.write("variant\tchrom\thap\trank\tcov_bp\n")
    for variant, idx in (("chr_cov", 1), ("fuse200k", 0)):
        for chrom in chroms:
            sc = chr_cov[chrom]
            for h, r in sorted(fused_rank[chrom][idx].items(), key=lambda x: x[1]):
                f.write(f"{variant}\t{chrom}\t{h}\t{r}\t{sc.get(h, 0)}\n")

# ---- genome-wide ranking: reciprocal-rank fusion over the groups ------------------------------------------
gchroms = [c for c in (a.groups or chroms) if c in fused_rank]
miss = [c for c in (a.groups or []) if c not in fused_rank]
if miss: print(f"[select] warning: no hits for {' '.join(miss)}", file=sys.stderr)
gscore = defaultdict(float)
for c in gchroms:
    for h, r in fused_rank[c][0].items(): gscore[h] += 1.0 / (K + r)
genome_order = sorted(gscore, key=lambda h: -gscore[h])
with open(a.out + ".genome_ranking.tsv", "w") as f:
    f.write("rank\thap\tscore\n")
    for i, h in enumerate(genome_order): f.write(f"{i + 1}\t{h}\t{gscore[h]:.6f}\n")

# ---- panels -------------------------------------------------------------------------------------------
for N in a.n:
    k = N - len(REFS)
    if k < 0: sys.exit(f"error: N={N} is smaller than the {len(REFS)} retained haplotypes")
    if a.mode in ("chr", "both"):
        with open(f"{a.out}.sift_chr_N{N}.tsv", "w") as f:
            for c in chroms:
                top = [h for h, _ in sorted(fused_rank[c][0].items(), key=lambda x: x[1])][:k]
                f.write(f"{c}\t{len(REFS) + len(top)}\t{','.join(REFS + top)}\n")
    if a.mode in ("genome", "both"):
        with open(f"{a.out}.sift_genome_N{N}.txt", "w") as f: f.write("\n".join(REFS + genome_order[:k]) + "\n")
print(f"[select] {len(panel)} haplotypes ranked, {len(chroms)} groups, {len(REFS)} retained; wrote {a.out}.rankings.tsv, "
      f"{a.out}.genome_ranking.tsv" + (f" and panels for N = {' '.join(map(str, a.n))} ({a.mode})" if a.n else ""),
      file=sys.stderr)
