#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fasta_toolkit.py -- FASTA utility functions for influenza multi-segment files.

Usage:
    python fasta_toolkit.py <command> [options]

Commands:
    unix            Strip Windows carriage returns and spaces after separator
    remove-doubles  Remove duplicate headers (keep first occurrence)
    fix-nuc         Replace non-nucleotide characters with 'n'
    split           Split into per-segment files
    concatenate     Concatenate segments per sample in canonical order (PB2 PB1 PA HA NP NA MP NS)
    subset-headers  Keep only selected fields from headers
    subsample       Filter by date, country, max-per-group, and/or similarity threshold

Run a command with --help for its options, e.g.:
    python fasta_toolkit.py subsample --help
"""

import argparse
import re
import sys
import os
import datetime
import random
from itertools import combinations

SEG_ORDER = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
NON_NUC   = re.compile(r"[^ACGTURYSWKMBDHVNacgturyswkmbdhvn?-]")


# -- I/O helpers ---------------------------------------------------------------

def read_fasta(path):
    """Yield (header_without_gt, sequence) tuples."""
    header, seq_lines = None, []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def write_fasta(records, path):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")


def auto_out(input_path, suffix):
    return os.path.splitext(input_path)[0] + suffix


# -- Functions -----------------------------------------------------------------

def to_unix(input_file, output_file, sep="|"):
    """
    Strip Windows carriage returns (\\r) and remove spaces directly after the separator.
    Reads and rewrites the file; operates line by line so large files are fine.
    """
    with open(input_file, "r", encoding="utf-8", errors="replace") as fin, \
         open(output_file, "w") as fout:
        for line in fin:
            line = line.replace("\r", "")
            line = re.sub(re.escape(sep) + r" +", sep, line)
            fout.write(line)
    print(f"Unix-converted -> {output_file}")


def remove_doubles(input_file, output_file):
    """
    Remove sequences with a duplicate header; keep the first occurrence.
    Comparison is on the full header string (excluding the leading '>').
    """
    seen = set()
    kept = removed = 0
    results = []
    for header, seq in read_fasta(input_file):
        if header in seen:
            removed += 1
        else:
            seen.add(header)
            results.append((header, seq))
            kept += 1
    write_fasta(results, output_file)
    print(f"Removed {removed} duplicate(s), kept {kept} -> {output_file}")


def fix_nucleotides(input_file, output_file):
    """Replace every character that is not a valid (degenerate) nucleotide with 'n'."""
    results = [(h, NON_NUC.sub("n", s)) for h, s in read_fasta(input_file)]
    write_fasta(results, output_file)
    print(f"Non-nucleotides replaced -> {output_file}")


def split_by_segment(input_file, output_prefix, sep="|", seg_index=-1):
    """
    Split a multi-segment FASTA into one file per segment.

    sep        : field separator in the header (default '|')
    seg_index  : 0-based position of the segment field; -1 = last field (default)
                 The segment name is taken as the part after the final '/' in that field.
    """
    base = output_prefix or os.path.splitext(input_file)[0]
    buckets = {}
    for header, seq in read_fasta(input_file):
        parts = header.split(sep)
        try:
            field = parts[seg_index]
        except IndexError:
            field = "unknown"
        segment = (field.split("/")[-1].split("_")[-1]).strip() or "unknown"
        buckets.setdefault(segment, []).append((header, seq))

    for segment, records in buckets.items():
        out_path = f"{base}_{segment}.fasta"
        write_fasta(records, out_path)
        print(f"  {segment}: {len(records)} sequences -> {out_path}")


def concatenate(input_file, output_file, sep="|", seg_index=-1):
    """
    Concatenate segments per sample in canonical influenza order:
    PB2 > PB1 > PA > HA > NP > NA > MP > NS.

    The segment field (seg_index) is stripped from the header to form the sample key.
    Segments missing for a sample are silently skipped.

    sep        : field separator in the header (default '|')
    seg_index  : 0-based position of the segment field; -1 = last field (default)
    """
    # group_key -> {segment_name: sequence}
    # group_key -> representative output header (fields minus segment field)
    sample_segs = {}
    sample_header = {}

    for header, seq in read_fasta(input_file):
        parts = header.split(sep)
        idx = seg_index if seg_index >= 0 else len(parts) + seg_index
        try:
            raw = parts[idx].strip()
            segment = raw.split("/")[-1].split("_")[-1]
        except IndexError:
            segment = "unknown"

        sample_parts = [p for i, p in enumerate(parts) if i != idx]

        # If the first remaining field looks like a bare accession (no "/" in
        # it) and there are at least 2 fields left, drop it from the group key
        # so that sequences from the same strain with different accessions
        # still merge correctly (common in GISAID multi-segment exports).
        if len(sample_parts) >= 2 and "/" not in sample_parts[0]:
            group_key = sep.join(sample_parts[1:])
        else:
            group_key = sep.join(sample_parts)

        sample_segs.setdefault(group_key, {})[segment] = seq
        if group_key not in sample_header:
            sample_header[group_key] = sep.join(sample_parts)

    base = os.path.splitext(os.path.basename(input_file))[0]
    out = output_file or "concatenated_{0}.fa".format(base)
    # Collect all segment names actually found in the file
    all_found_segs = set()
    for segs in sample_segs.values():
        all_found_segs.update(segs.keys())
    known = [s for s in SEG_ORDER if s in all_found_segs]
    unknown = sorted(all_found_segs - set(SEG_ORDER))
    if unknown:
        print("WARNING: segment name(s) not in canonical order, appended at end: {0}".format(unknown))
    order = known + unknown

    written = skipped = 0
    with open(out, "w") as fh:
        for group_key, segs in sample_segs.items():
            concat_seq = "".join(segs[s] for s in order if s in segs)
            if not concat_seq:
                skipped += 1
                continue
            fh.write(">{0}\n{1}\n".format(sample_header[group_key], concat_seq))
            written += 1
    if skipped:
        print("WARNING: {0} sample(s) produced empty sequences and were skipped.".format(skipped))
    print("Concatenated {0} samples -> {1}".format(written, out))


def subset_headers(input_file, output_file, indices, sep="|"):
    """
    Keep only the fields at the given indices from pipe-delimited headers.

    indices : list of 0-based integers, e.g. [0, 2]
    sep     : field separator (default '|')
    """
    results = []
    for header, seq in read_fasta(input_file):
        parts = header.split(sep)
        selected = [parts[i] for i in indices if i < len(parts)]
        if not selected:
            continue
        new_header = sep.join(selected)
        # Ensure the '>' is preserved correctly
        if not new_header.startswith(">"):
            pass  # read_fasta already strips '>'
        results.append((new_header, seq))
    write_fasta(results, output_file)
    print(f"Header subset -> {output_file}")


def subsample(
    input_file,
    output_file,
    sep="|",
    # Header field positions
    sample_field=1,       # field containing "type/country/strain" slash-info
    date_field=2,         # field containing the date (YYYY-MM-DD or YYYY)
    country_slash_pos=2,  # position within the slash-split sample_field
    # Date filtering
    date_after=None,      # keep sequences with date >= this (YYYY-MM-DD)
    date_before=None,     # keep sequences with date <= this (YYYY-MM-DD)
    # Country filtering
    include_countries=None,   # whitelist: only keep these countries
    exclude_countries=None,   # blacklist: always remove these countries
    keep_countries=None,      # exempt from similarity removal (keep all)
    # Per-group cap (applied after date filter, before similarity)
    max_per_group=None,   # max sequences per (country, year, month) group
    # Similarity filtering
    similarity_threshold=None,  # float 0-100; remove near-identical seqs within group
):
    """
    Subsample a FASTA in up to four sequential steps:

    1. Date filter   -- drop sequences outside [date_after, date_before]
    2. Country filter -- apply include/exclude country lists
    3. Per-group cap  -- randomly keep at most max_per_group sequences per
                        (country, year, month) group
    4. Similarity     -- within each group, remove sequences >= similarity_threshold %
                        similar to another (longer sequence wins); countries listed in
                        keep_countries are exempt from this step

    Header format assumed (pipe-separated, 0-based):
        field[0]            : accession / ID
        field[sample_field] : type/country/strain  (slash-separated;
                              country at slash position country_slash_pos)
        field[date_field]   : YYYY-MM-DD  (or YYYY)

    All filtering steps are optional -- only those with a non-None argument run.
    """
    try:
        from Bio import Align
        have_biopython = True
    except ImportError:
        have_biopython = False
        if similarity_threshold is not None:
            sys.exit("Biopython is required for similarity filtering. "
                     "Install with: pip install biopython")

    # -- helpers --------------------------------------------------------------

    def parse_date(s):
        """Parse YYYY-MM-DD or YYYY into a date object; return None on failure."""
        s = s.strip()
        for fmt in ("%Y-%m-%d", "%Y-%m", "%Y"):
            try:
                return datetime.datetime.strptime(s, fmt).date()
            except ValueError:
                continue
        return None

    def extract_meta(header):
        """Return (country, date_obj) from a header string."""
        parts = header.split(sep)
        # Country
        try:
            slash_parts = parts[sample_field].split("/")
            country = (slash_parts[country_slash_pos]
                       if len(slash_parts) > country_slash_pos
                       else parts[sample_field])
        except IndexError:
            country = header
        # Date
        try:
            date_obj = parse_date(parts[date_field])
        except IndexError:
            date_obj = None
        return country.strip(), date_obj

    def group_key(country, date_obj):
        if date_obj:
            return (country, str(date_obj.year), f"{date_obj.month:02d}")
        return (country, "unknown", "unknown")

    def calc_similarity(seq1, seq2):
        aligner = Align.PairwiseAligner()
        score = aligner.align(seq1, seq2).score
        return score / max(len(seq1), len(seq2)) * 100

    # -- load -----------------------------------------------------------------

    print(f"[subsample] Reading {input_file} ...")
    records = list(read_fasta(input_file))
    print(f"  {len(records)} sequences loaded  ({datetime.datetime.now():%H:%M:%S})")

    d_after  = parse_date(date_after)  if date_after  else None
    d_before = parse_date(date_before) if date_before else None
    inc_set  = {c.strip() for c in include_countries} if include_countries else None
    exc_set  = {c.strip() for c in exclude_countries} if exclude_countries else set()
    keep_set = {c.strip() for c in keep_countries}    if keep_countries    else set()

    # -- step 1 : date filter -------------------------------------------------

    if d_after or d_before:
        before = len(records)
        def date_ok(header):
            _, d = extract_meta(header)
            if d is None:
                return True          # no date -> keep (can't judge)
            if d_after  and d < d_after:  return False
            if d_before and d > d_before: return False
            return True
        records = [(h, s) for h, s in records if date_ok(h)]
        print(f"  After date filter:    {len(records)}  (removed {before - len(records)})")

    # -- step 2 : country filter ----------------------------------------------

    if inc_set or exc_set:
        before = len(records)
        def country_ok(header):
            c, _ = extract_meta(header)
            if exc_set and c in exc_set: return False
            if inc_set and c not in inc_set: return False
            return True
        records = [(h, s) for h, s in records if country_ok(h)]
        print(f"  After country filter: {len(records)}  (removed {before - len(records)})")

    # -- step 3 : per-group cap -----------------------------------------------

    if max_per_group is not None:
        before = len(records)
        groups = {}
        for h, s in records:
            c, d = extract_meta(h)
            k = group_key(c, d)
            groups.setdefault(k, []).append((h, s))
        records = []
        for k, grp in groups.items():
            if len(grp) > max_per_group:
                grp = random.sample(grp, max_per_group)
            records.extend(grp)
        print(f"  After per-group cap ({max_per_group}): {len(records)}  (removed {before - len(records)})")

    # -- step 4 : similarity filter -------------------------------------------

    if similarity_threshold is not None:
        before = len(records)
        # Build groups
        groups = {}
        for i, (h, s) in enumerate(records):
            c, d = extract_meta(h)
            k = group_key(c, d)
            groups.setdefault(k, []).append(i)

        remove_idx = set()
        n_removed  = 0
        for (country, year, month), indices in groups.items():
            if country in keep_set:
                continue
            for i, j in combinations(indices, 2):
                if i in remove_idx or j in remove_idx:
                    continue
                sim = calc_similarity(records[i][1], records[j][1])
                if sim >= similarity_threshold:
                    # keep the longer sequence
                    drop = j if len(records[i][1]) >= len(records[j][1]) else i
                    remove_idx.add(drop)
                    n_removed += 1
                    if n_removed % 500 == 0:
                        print(f"    similarity: removed {n_removed} so far ...")

        records = [(h, s) for idx, (h, s) in enumerate(records) if idx not in remove_idx]
        pct = n_removed / before * 100 if before else 0
        print(f"  After similarity filter (>={similarity_threshold}%): "
              f"{len(records)}  (removed {n_removed} = {pct:.1f}%)")

    # -- write -----------------------------------------------------------------

    out = output_file or auto_out(input_file, "_subsampled.fasta")
    write_fasta(records, out)
    print(f"  Done -> {out}  ({datetime.datetime.now():%H:%M:%S})")


# -- CLI -----------------------------------------------------------------------

def _sep_arg(p):
    p.add_argument("--sep", default="|",
                   help="Field separator in headers (default '|')")

def _seg_arg(p):
    p.add_argument("--seg-index", dest="seg_index", type=int, default=-1,
                   help="0-based index of the segment field (default: -1 = last field)")

def make_parser():
    parser = argparse.ArgumentParser(
        prog="fasta_toolkit.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", metavar="command")
    sub.required = True

    # unix
    p = sub.add_parser("unix", help="Strip \\r and spaces after separator")
    p.add_argument("input");  p.add_argument("-o", "--output", default=None)
    _sep_arg(p)

    # remove-doubles
    p = sub.add_parser("remove-doubles", help="Remove duplicate headers")
    p.add_argument("input");  p.add_argument("-o", "--output", default=None)

    # fix-nuc
    p = sub.add_parser("fix-nuc", help="Replace non-nucleotide characters with 'n'")
    p.add_argument("input");  p.add_argument("-o", "--output", default=None)

    # split
    p = sub.add_parser("split", help="Split into one file per segment")
    p.add_argument("input")
    p.add_argument("--prefix", default=None, help="Output filename prefix (default: input basename)")
    _sep_arg(p); _seg_arg(p)

    # concatenate
    p = sub.add_parser("concatenate", help="Concatenate segments per sample in canonical order")
    p.add_argument("input");  p.add_argument("-o", "--output", default=None)
    _sep_arg(p); _seg_arg(p)

    # subset-headers
    p = sub.add_parser("subset-headers", help="Keep only selected header fields")
    p.add_argument("input");  p.add_argument("output")
    p.add_argument("indices", nargs="+", type=int, help="0-based field indices to keep")
    _sep_arg(p)

    # subsample
    p = sub.add_parser(
        "subsample",
        help="Filter by date, country, per-group cap, and/or similarity",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Subsample a FASTA through up to four optional, sequential steps:

  1. Date filter       --after / --before
  2. Country filter    --include-countries / --exclude-countries
  3. Per-group cap     --max-per-group  (random; reproducible with --seed)
  4. Similarity        --similarity  (within country/year/month groups)

Header format expected (pipe-separated, 0-based fields):
  field[--sample-field]  type/country/strain  (country at --country-pos in slash-split)
  field[--date-field]    YYYY-MM-DD

Example -- keep only 2020-onwards European sequences, max 10/group, >99% similar removed,
but always keep all Belgian sequences:
  python fasta_toolkit.py subsample in.fasta -o out.fasta \\
      --after 2020-01-01 \\
      --include-countries Netherlands Belgium Germany France \\
      --max-per-group 10 \\
      --similarity 99 \\
      --keep-countries Belgium
""",
    )
    p.add_argument("input")
    p.add_argument("-o", "--output", default=None)
    _sep_arg(p)
    # Header layout
    p.add_argument("--sample-field",  dest="sample_field",      type=int, default=1,
                   help="Field index containing type/country/strain (default: 1)")
    p.add_argument("--date-field",    dest="date_field",        type=int, default=2,
                   help="Field index containing the date (default: 2)")
    p.add_argument("--country-pos",   dest="country_slash_pos", type=int, default=2,
                   help="Position within slash-split sample field where country sits (default: 2)")
    # Date
    p.add_argument("--after",  dest="date_after",  default=None, metavar="YYYY-MM-DD",
                   help="Keep sequences with date >= this")
    p.add_argument("--before", dest="date_before", default=None, metavar="YYYY-MM-DD",
                   help="Keep sequences with date <= this")
    # Country
    p.add_argument("--include-countries", dest="include_countries", nargs="+", default=None,
                   metavar="COUNTRY",
                   help="Whitelist: only keep sequences from these countries")
    p.add_argument("--exclude-countries", dest="exclude_countries", nargs="+", default=None,
                   metavar="COUNTRY",
                   help="Blacklist: always remove sequences from these countries")
    p.add_argument("--keep-countries",    dest="keep_countries",    nargs="+", default=None,
                   metavar="COUNTRY",
                   help="Exempt these countries from similarity removal (keep all their sequences)")
    # Per-group cap
    p.add_argument("--max-per-group", dest="max_per_group", type=int, default=None,
                   metavar="N",
                   help="Randomly keep at most N sequences per (country, year, month) group")
    p.add_argument("--seed", dest="seed", type=int, default=None,
                   help="Random seed for reproducible per-group sampling")
    # Similarity
    p.add_argument("--similarity", dest="similarity_threshold", type=float, default=None,
                   metavar="THRESHOLD",
                   help="Remove sequences >= THRESHOLD %% similar within each group (e.g. 99)")

    return parser


def main():
    parser = make_parser()
    args = parser.parse_args()

    if args.command == "unix":
        out = args.output or auto_out(args.input, "_unix.fasta")
        to_unix(args.input, out, sep=args.sep)

    elif args.command == "remove-doubles":
        out = args.output or auto_out(args.input, "_nodups.fasta")
        remove_doubles(args.input, out)

    elif args.command == "fix-nuc":
        out = args.output or auto_out(args.input, "_fixnuc.fasta")
        fix_nucleotides(args.input, out)

    elif args.command == "split":
        split_by_segment(args.input, args.prefix, sep=args.sep, seg_index=args.seg_index)

    elif args.command == "concatenate":
        concatenate(args.input, args.output, sep=args.sep, seg_index=args.seg_index)

    elif args.command == "subset-headers":
        subset_headers(args.input, args.output, args.indices, sep=args.sep)

    elif args.command == "subsample":
        if args.seed is not None:
            random.seed(args.seed)
        subsample(
            input_file=args.input,
            output_file=args.output,
            sep=args.sep,
            sample_field=args.sample_field,
            date_field=args.date_field,
            country_slash_pos=args.country_slash_pos,
            date_after=args.date_after,
            date_before=args.date_before,
            include_countries=args.include_countries,
            exclude_countries=args.exclude_countries,
            keep_countries=args.keep_countries,
            max_per_group=args.max_per_group,
            similarity_threshold=args.similarity_threshold,
        )


if __name__ == "__main__":
    main()
