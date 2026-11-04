# SIFT

SIFT indexes a pangenome parsed with [PFSP](https://github.com/<user>/PFSP) and finds pseudo-MEMs between a query
and the pangenome's haplotypes. From the pseudo-MEMs of an assembly's contigs, `scripts/select_panels.py` ranks the
haplotypes closest to the assembly, without population labels, and selects a small donor panel for assembly
patching (for example with ImpuT2T).

## Build

Requirements: a C++20 compiler, CMake ≥ 3.20, zlib, OpenMP, and git. sdsl-lite, CLI11, kseq and MurmurHash3 are
downloaded at configure time at pinned versions; PFSP is included as a git submodule (`external/PFSP`).

```bash
git clone --recursive https://github.com/<user>/SIFT.git
cd SIFT
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j        # -> build/sift
```

To use a PFSP checkout elsewhere, configure with `-DPFSP_DIR=/path/to/PFSP`.

## Usage

Haplotypes are tagged by their PanSN name (`sample#haplotype#contig`). The pangenome can be indexed whole or in
parts; one index per chromosome keeps the memory of the build down. The examples use chromosome 1.

**1. Parse** the pangenome with PFSP:
```bash
PFSP chr1.pangenome.fa -o panel/chr1 -k 200 -s 20 --seed 0
```

**2. Build** the index: forward and reverse run-length compressed suffix arrays of the parse, and the haplotype
tag array (`--tag-sample s` stores every s-th tag and recovers the rest with at most s − 1 Psi steps). Writes
`index/chr1.fwd.rlcsa`, `index/chr1.rev.rlcsa` and `index/chr1.tags`. Construction runs in memory.
```bash
build/sift build panel/chr1 -o index/chr1 --tag-by pansn-hap --tag-sample 8
```

**3. Query** a FASTA/FASTQ file (optionally gzipped) of sequences, for example the assembly contigs placed on
chromosome 1:
```bash
build/sift query contigs/chr1.fa -d panel/chr1 -i index/chr1 -L 1000 -t 8 -o /dev/stdout |
    grep '^>' | gzip > hits/chr1.L1000.hits.gz
```
Output: one FASTA record per pseudo-MEM, `>name:start-end tags=sample#hap,...`. With `-t` > 1, records come out in
the order the threads finish.

`sift query` searches the forward strand of each sequence only. To match both strands of a contig, write its
reverse complement as a second record. Naming the two records `assembly_hap|contig|+` and `assembly_hap|contig|-`
lets `select_panels.py` score them as one contig. For example, with [seqkit](https://bioinf.shenwei.me/seqkit/):
```bash
seqkit replace -p '^(\S+).*' -r 'hap1|$1|+' asm.hap1.chr1.fa  > contigs/chr1.fa
seqkit seq -t dna -r -p asm.hap1.chr1.fa | seqkit replace -p '^(\S+).*' -r 'hap1|$1|-' >> contigs/chr1.fa
```

**4. Select panels** from the hits:
```bash
python3 scripts/select_panels.py hits/chr*.L1000.hits.gz --contigs contigs/ \
    --exclude HG002 --always CHM13 GRCh38 --n 10 25 50 --mode both --out HG002.select
```
- Each hits file is ranked as its own group. A file whose name contains a chromosome (`chr1` … `chrY`) belongs to
  that chromosome; any other file is a group named after the file. A single whole-genome hits file is one group.
- `--exclude` drops samples or haplotypes from the ranking (the query sample, its relatives). `--always` keeps
  samples or haplotypes, such as the references, in every panel without ranking them. Names can be samples
  (`HG002`, all of its haplotypes), haplotypes (`HG002.1`) or PanSN (`HG002#1`). `--exclude-file` and
  `--always-file` read the same names from a file, one per line.
- `--contigs` (FASTA files, or folders holding `<group>.fa`) gives exact sequence lengths. Without it, lengths are
  estimated from the hits.

Each panel holds the `--always` haplotypes and the top-ranked haplotypes, N in total:
- `--mode chr`: a panel per group (SIFT-CHR), in `HG002.select.sift_chr_N<N>.tsv`.
- `--mode genome`: one panel for the whole genome (SIFT-Genome), in `HG002.select.sift_genome_N<N>.txt`.

The rankings are in `HG002.select.rankings.tsv` (per group) and `HG002.select.genome_ranking.tsv`.
`python3 scripts/select_panels.py --help` describes the scoring.

## Acknowledgements

The run-length compressed suffix arrays build on the EFRank structure of Travis Gagie's
[KATKA](https://github.com/TravisGagie/KATKA).

## Citation

## License

GNU General Public License v3.0 (see `LICENSE`).
