# Influenza FASTA Toolkit

A collection of Python scripts for cleaning, manipulating, and subsampling influenza multi-segment FASTA files.

---

## Scripts

| Script | Purpose |
|---|---|
| `fasta_toolkit.py` | All-in-one toolkit: clean, split, concatenate, subsample, and more |
| `trimmer.py` | Trim sequences to coding region (first ATG → last stop codon) |
| `welke_mis_ik.py` | Find FASTA headers that are absent from a reference Excel file |

---

## Requirements

```bash
pip install biopython pandas openpyxl
```

`biopython` is only required for the `subsample --similarity` filter in `fasta_toolkit.py`, and for both `trimmer.py` and `welke_mis_ik.py`. The other toolkit commands run on the Python standard library alone.

---

## Expected header format

Scripts assume pipe-separated headers in INSDC/GISAID style (0-based field indices):

```
>0_accession|1_type/country/strain|2_YYYY-MM-DD|3_segment
```

Example:
```
>123456|A/Belgium/UZ001/2023|2023-11-14|PB2
```

Where the segment field is the last field by default. All field positions are configurable — see `--sep`, `--seg-index`, `--sample-field`, `--date-field`, and `--country-pos`.

---

## fasta_toolkit.py

```
python fasta_toolkit.py <command> [options]
```

Run any command with `--help` to see its full option list:
```bash
python fasta_toolkit.py subsample --help
```

---

### unix

Strip Windows carriage returns (`\r`) and remove spaces that appear directly after the field separator.

```bash
python fasta_toolkit.py unix <input> [-o output] [--sep SEP]
```

| Argument | Default | Description |
|---|---|---|
| `input` | — | Input FASTA file |
| `-o` / `--output` | `<input>_unix.fasta` | Output file |
| `--sep` | `\|` | Field separator to clean spaces after |

```bash
python fasta_toolkit.py unix sequences.fasta
python fasta_toolkit.py unix sequences.fasta -o sequences_clean.fasta
```

---

### remove-doubles

Remove sequences with a duplicate header, keeping the first occurrence. Comparison is on the full header string.

```bash
python fasta_toolkit.py remove-doubles <input> [-o output]
```

| Argument | Default | Description |
|---|---|---|
| `input` | — | Input FASTA file |
| `-o` / `--output` | `<input>_nodups.fasta` | Output file |

```bash
python fasta_toolkit.py remove-doubles sequences.fasta
```

---

### fix-nuc

Replace any character that is not a valid IUPAC nucleotide (including degenerate bases) with `n`.

Valid characters retained: `A C G T U R Y S W K M B D H V N ? -` (upper and lower case).

```bash
python fasta_toolkit.py fix-nuc <input> [-o output]
```

| Argument | Default | Description |
|---|---|---|
| `input` | — | Input FASTA file |
| `-o` / `--output` | `<input>_fixnuc.fasta` | Output file |

```bash
python fasta_toolkit.py fix-nuc sequences.fasta
```

---

### split

Split a multi-segment FASTA into one file per segment. The segment name is read from a configurable field in the header.

```bash
python fasta_toolkit.py split <input> [--prefix PREFIX] [--sep SEP] [--seg-index N]
```

| Argument | Default | Description |
|---|---|---|
| `input` | — | Input FASTA file |
| `--prefix` | Input basename | Prefix for output filenames |
| `--sep` | `\|` | Header field separator |
| `--seg-index` | `-1` (last field) | 0-based index of the segment field |

Output files are named `<prefix>_<segment>.fasta`, e.g. `sequences_PB2.fasta`.

```bash
# Standard pipe-separated headers, segment in last field
python fasta_toolkit.py split sequences.fasta

# Segment in field index 2, slash-separated headers
python fasta_toolkit.py split sequences.fasta --sep / --seg-index 2

# Custom output prefix
python fasta_toolkit.py split sequences.fasta --prefix results/2024
```

---

### concatenate

Concatenate all segments per sample into a single sequence in canonical influenza order: **PB2 › PB1 › PA › HA › NP › NA › MP › NS**.

Samples are grouped by their header with the segment field removed. Segments missing for a sample are silently skipped.

```bash
python fasta_toolkit.py concatenate <input> [-o output] [--sep SEP] [--seg-index N]
```

| Argument | Default | Description |
|---|---|---|
| `input` | — | Input FASTA file containing all segments |
| `-o` / `--output` | `concatenated_<input>.fa` | Output file |
| `--sep` | `\|` | Header field separator |
| `--seg-index` | `-1` (last field) | 0-based index of the segment field |

```bash
python fasta_toolkit.py concatenate all_segments.fasta -o full_genomes.fa
```

---

### subset-headers

Keep only the specified fields from each header, discarding the rest. Fields are identified by their 0-based index in the separator-split header.

```bash
python fasta_toolkit.py subset-headers <input> <output> <index> [index ...] [--sep SEP]
```

| Argument | Default | Description |
|---|---|---|
| `input` | — | Input FASTA file |
| `output` | — | Output file (required) |
| `indices` | — | One or more 0-based field indices to keep |
| `--sep` | `\|` | Header field separator |

```bash
# Keep accession (field 0) and date (field 2) only
python fasta_toolkit.py subset-headers sequences.fasta out.fasta 0 2

# Keep fields 0, 1, and 3
python fasta_toolkit.py subset-headers sequences.fasta out.fasta 0 1 3
```

---

### subsample

Filter a FASTA through up to four optional, sequential steps. All filters are independent — use any combination.

```bash
python fasta_toolkit.py subsample <input> [-o output] [--sep SEP]
    [--sample-field N] [--date-field N] [--country-pos N]
    [--after YYYY-MM-DD] [--before YYYY-MM-DD]
    [--include-countries C [C ...]] [--exclude-countries C [C ...]]
    [--keep-countries C [C ...]]
    [--max-per-group N] [--seed N]
    [--similarity THRESHOLD]
```

#### Steps (applied in this order)

**Step 1 — Date filter**

Keep only sequences whose date falls within the specified range. Sequences with unparseable or missing dates are kept.

| Argument | Description |
|---|---|
| `--after YYYY-MM-DD` | Keep sequences with date ≥ this |
| `--before YYYY-MM-DD` | Keep sequences with date ≤ this |

**Step 2 — Country filter**

| Argument | Description |
|---|---|
| `--include-countries` | Whitelist: only keep sequences from these countries |
| `--exclude-countries` | Blacklist: remove sequences from these countries |

Both can be used together. Exclusion is applied first.

**Step 3 — Per-group cap**

Randomly subsample within each (country, year, month) group down to at most N sequences. Applied before similarity filtering so the similarity step works on a smaller set.

| Argument | Description |
|---|---|
| `--max-per-group N` | Maximum sequences per (country, year, month) group |
| `--seed N` | Random seed for reproducible results |

**Step 4 — Similarity filter**

Within each (country, year, month) group, compare all pairs of sequences using pairwise alignment. When two sequences exceed the threshold, the shorter one is removed. Requires Biopython.

| Argument | Description |
|---|---|
| `--similarity THRESHOLD` | Remove sequences ≥ THRESHOLD % similar within each group (e.g. `99`) |
| `--keep-countries C [C ...]` | Exempt these countries from similarity removal entirely |

#### Header layout options

| Argument | Default | Description |
|---|---|---|
| `--sep` | `\|` | Field separator |
| `--sample-field N` | `1` | Field index containing `type/country/strain` |
| `--date-field N` | `2` | Field index containing the date |
| `--country-pos N` | `2` | Position of country within the slash-split sample field |

#### Examples

```bash
# Keep only sequences from 2020 onwards
python fasta_toolkit.py subsample sequences.fasta -o out.fasta --after 2020-01-01

# Keep only Belgian and Dutch sequences
python fasta_toolkit.py subsample sequences.fasta -o out.fasta \
    --include-countries Belgium Netherlands

# Remove sequences that are ≥99% similar within each country/year/month group,
# but always keep all Belgian sequences
python fasta_toolkit.py subsample sequences.fasta -o out.fasta \
    --similarity 99 --keep-countries Belgium

# Full pipeline: date range + country whitelist + max 10/group + similarity filter
python fasta_toolkit.py subsample sequences.fasta -o out.fasta \
    --after 2020-01-01 --before 2024-12-31 \
    --include-countries Belgium Netherlands Germany France \
    --max-per-group 10 --seed 42 \
    --similarity 99 --keep-countries Belgium
```

---

## trimmer.py

Trim each sequence in a FASTA to its coding region by finding the first start codon (`ATG`) and the last stop codon (`TAA`, `TAG`, or `TGA`). Sequences where either codon is not found are written untrimmed.

**Requires:** `biopython`

```bash
python trimmer.py <input.fasta> <output.fasta>
```

| Argument | Description |
|---|---|
| `input.fasta` | Input FASTA file |
| `output.fasta` | Output FASTA file with trimmed sequences |

```bash
python trimmer.py sequences.fasta sequences_trimmed.fasta
```

> **Note:** The script searches for the last occurrence of any stop codon across the whole sequence, which works well for full coding sequences but may give unexpected results if your sequences contain internal stop codons (e.g. pseudogenes or sequencing errors).

---

## welke_mis_ik.py

Compare the headers in a FASTA file against the first column of an Excel file and print any headers that are missing from the Excel. Useful for checking whether all sequences in a dataset are represented in a metadata sheet.

**Requires:** `biopython`, `pandas`, `openpyxl`

Before comparison, spaces in the Excel column are replaced with underscores and any leading `>` is stripped, to match FASTA header formatting.

```bash
python welke_mis_ik.py <input.fasta> <reference.xlsx>
```

| Argument | Description |
|---|---|
| `input.fasta` | FASTA file whose headers you want to check |
| `reference.xlsx` | Excel file with reference identifiers in the first column |

```bash
python welke_mis_ik.py sequences.fasta metadata.xlsx
```

Output is printed to stdout:
```
Headers from the multi-FASTA file not found in the first column of the Excel file:
A/Belgium/UZ042/2023
A/Netherlands/NL099/2022
```

---

## Typical workflow

```bash
# 1. Fix line endings and clean up headers
python fasta_toolkit.py unix raw.fasta -o step1.fasta

# 2. Remove duplicate sequences
python fasta_toolkit.py remove-doubles step1.fasta -o step2.fasta

# 3. Replace non-nucleotide characters
python fasta_toolkit.py fix-nuc step2.fasta -o step3.fasta

# 4. Subsample: keep 2020–2024, max 20 per country/month, remove near-identical
python fasta_toolkit.py subsample step3.fasta -o step4.fasta \
    --after 2020-01-01 --before 2024-12-31 \
    --max-per-group 20 --seed 1 \
    --similarity 99 --keep-countries Belgium

# 5. Concatenate all segments into full genomes
python fasta_toolkit.py concatenate step4.fasta -o full_genomes.fa

# 6. Check all sequences are in your metadata sheet
python welke_mis_ik.py full_genomes.fa metadata.xlsx
```
