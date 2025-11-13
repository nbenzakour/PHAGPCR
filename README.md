# PHAGPCR - qPCR Primer Design for Phage Detection

Design and evaluate specific qPCR primers for phage detection in patient serum with comprehensive specificity screening.

![Workflow](data/Figma_2023-03-07.png)

---

## Abstract

Screening phages directly from a patient's serum using qPCR methods can be incredibly valuable in identifying phage infections and monitoring their response to treatment. To ensure the accuracy and specificity of qPCR results, it is important to carefully design primers that can detect the target phages with minimal cross-reactivity with genetic material that can be present in a serum sample, namely bacterial, viral and human DNA.

Here we present **PHAGPCR**, a pipeline built around Primer3 and MFEprimer, that can help identify suitable primer sequences and assess their specificity against local and public databases representing the phageome, microbiome and human genome. By combining these, researchers and clinicians can minimize the risk of false positive results and obtain a clearer picture of the phage populations present in a patient's serum.

Furthermore, this approach can be used to predict primers compatible with tracking multiple phages from within the same serum sample, provided that unique regions exist and suitable primers can be predicted. This can be particularly useful in cases where a combination of phages is administered to a patient for therapeutic treatment.

---

## Installation

### Requirements

- Python 3.9+
- MFEprimer 3.2 (local binary in `bin/mfeprimer`)
- BLAST+ (optional, for custom databases)

### Python Dependencies

Install required packages:

```bash
pip install -r requirements.txt
```

**Dependencies:**
- biopython >= 1.81
- pandas >= 2.0.0
- primer3-py >= 2.0.0
- termcolor >= 2.0.0

### MFEprimer Setup

Place the MFEprimer executable in `bin/mfeprimer` or update the `MFEPRIMER_BIN` constant in `phagpcr.py`.

---

## Usage

### Basic Commands

#### 1. Single Target Primer Design (Runtype 1)

Design primers for a specific gene or sequence:

```bash
python phagpcr.py \
  -f data/targets/CDS_0033_primase.fna \
  -o results/target_primers \
  -r 1 \
  -nb 3
```

#### 2. Cocktail/Multiplex Detection (Runtype 2)

Design primers for multiple phages with compatibility checking:

```bash
python phagpcr.py \
  -f multi_phage_targets.fasta \
  -o results/cocktail_primers \
  -r 2 \
  -nb 5
```

#### 3. Unique Region Detection (Runtype 3)

**Status:** Planned - Not yet implemented
See `test/planning/Phase_3_planning.md` for implementation details.

### Advanced Options

#### With BLAST Database Screening

```bash
python phagpcr.py \
  -f input.fasta \
  -o results/ \
  -b data/blastdb/meta_ILL_ONT_20250521_combined.fna \
  -r 1 \
  -nb 3
```

#### With Human Genome Screening

```bash
python phagpcr.py \
  -f input.fasta \
  -o results/ \
  -hg \
  -r 1 \
  -nb 3
```

**Note:** Human genome screening requires hg38 database at default path or specify custom path in code.

---

## Command Line Arguments

| Argument | Short | Type | Default | Description |
|----------|-------|------|---------|-------------|
| `--fasta_file` | `-f` | str | Required | Path to input FASTA file |
| `--output_dir` | `-o` | str | Required | Output directory for results |
| `--blast_db` | `-b` | str | `data/blastdb/...` | BLAST database for screening |
| `--human_genome` | `-hg` | flag | False | Screen against human genome (hg38) |
| `--sequences_dir` | `-s` | str | None | Directory with sequences for screening |
| `--runtype` | `-r` | int | 1 | Run type: 1=target, 2=cocktail, 3=unique |
| `--primer_nb` | `-nb` | int | 1 | Number of primer pairs per target |
| `--tm_optimal` | `-tm` | int | 60 | Optimal Tm for primers (°C) |
| `--size_optimal` | `-sz` | int | 20 | Optimal primer size (bp) |
| `--kit` | `-k` | str | "" | qPCR kit (Quantinova/Invitrogen) |
| `--tm_specificity` | `-tm2` | int | 40 | Min Tm for MFEprimer matches (°C) |

---

## Workflows

### Workflow 1: TARGET qPCR (runtype=1)

**Status:** ✅ Fully Functional

**Purpose:** Design primers for a specific target gene or sequence

**Steps:**
1. Input: Single target sequence (FASTA)
2. Primer3 design
3. MFEprimer screening against custom database
4. MFEprimer screening against human genome (optional)

**Output:**
- `prefiltered_primers_file.csv` - All designed primers
- `prefiltered_primers.fna` - Primers in FASTA format
- `mfe_screening_results_against_local_database_summary.csv` - Specificity results
- `mfe_spec_results_hg38_summary.csv` - Human genome hits (if `-hg`)

---

### Workflow 2: COCKTAIL Detection (runtype=2)

**Status:** ✅ Fully Functional

**Purpose:** Design compatible primers for multiple phages (multiplex qPCR)

**Steps:**
1. Input: Multiple target sequences (multiFASTA)
2. Primer3 design for each target
3. MFEprimer screening against custom database
4. MFEprimer screening against human genome (optional)
5. MFEprimer dimer analysis (inter-primer compatibility)

**Output:**
- Same as Workflow 1, plus:
- `mfe_dimer_results_prefiltered_primers_summary.csv` - Dimer analysis
  - "Problematic" column identifies inter-target dimers

---

### Workflow 3: UNIQUE TARGET (runtype=3)

**Status:** ❌ Not Implemented (Planned)

**Purpose:** Find unique regions in phage genome when specific genes lack uniqueness

See `test/planning/Phase_3_planning.md` for detailed implementation plan.

---

## Output Files

### Standard Output (All Runtypes)

| File | Description |
|------|-------------|
| `prefiltered_primers_file.csv` | All primers with Tm, GC%, penalties, etc. |
| `prefiltered_primers.fna` | Primers in FASTA format |
| `mfe_screening_results_against_local_database_summary.csv` | Specificity against custom DB |
| `mfe_spec_results_hg38_summary.csv` | Human genome screening results |

### Cocktail-Specific (Runtype 2)

| File | Description |
|------|-------------|
| `mfe_dimer_results_prefiltered_primers_summary.csv` | Primer-primer interactions |

### Additional Files

| Directory/File | Description |
|----------------|-------------|
| `indiv_results/` | Individual MFEprimer results per primer pair |

---

## Gene Extractor Tool

Extract specific genes from GenBank annotation files.

### Usage

**Single file:**
```bash
python bz_gene_extractor.py \
  -i data/genbank/phage.gbk \
  -w "primase" \
  -o extracted_genes/ \
  -sk 500
```

**Directory with fallback:**
```bash
python bz_gene_extractor.py \
  -d data/genbank/ \
  -w "DNA polymerase" \
  -w2 "polymerase" \
  -o extracted_genes/ \
  -sk 200
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-i, --input` | Input GenBank file |
| `-d, --directory` | Input directory (mutually exclusive with `-i`) |
| `-w, --word` | Primary keyword to search in product qualifiers |
| `-w2, --word2` | Fallback keyword if primary not found |
| `-o, --outdir` | Output directory for extracted sequences |
| `-sk, --skip` | Minimum sequence length (bp) to extract |

---

## Examples

### Example 1: Basic Primer Design

```bash
python phagpcr.py \
  -f data/targets/CDS_0033_primase.fna \
  -o test_output \
  -r 1 \
  -nb 3
```

**Result:** 3 primer pairs designed, no screening

### Example 2: Full Screening Pipeline

```bash
python phagpcr.py \
  -f data/targets/CDS_0033_primase.fna \
  -o test_output \
  -b data/blastdb/meta_ILL_ONT_20250521_combined.fna \
  -hg \
  -r 1 \
  -nb 5
```

**Result:** 5 primer pairs with full specificity screening

### Example 3: Multiplex Design

```bash
python phagpcr.py \
  -f cocktail_targets.fasta \
  -o cocktail_results \
  -b data/blastdb/meta_ILL_ONT_20250521_combined.fna \
  -r 2 \
  -nb 3
```

**Result:** Compatible primers for multiplex detection

---

## Troubleshooting

### Common Issues

#### 1. MFEprimer Not Found

**Error:** `FileNotFoundError: [Errno 2] No such file or directory: 'bin/mfeprimer'`

**Solution:** Ensure MFEprimer binary is in `bin/mfeprimer` or update `MFEPRIMER_BIN` constant

#### 2. No Primers Designed

**Error:** `ValueError: No primers designed for <sequence>`

**Possible Causes:**
- Sequence too short (<200bp)
- GC content too extreme
- Repetitive sequence

**Solutions:**
- Adjust Tm parameters (`-tm`)
- Adjust primer size (`-sz`)
- Try different target region

#### 3. BLAST Database Errors

**Error:** `subprocess.CalledProcessError: ... mfeprimer index ...`

**Solutions:**
- Ensure database file exists
- Check database format (FASTA)
- Verify file permissions

#### 4. Human Genome Path Not Found

**Error:** Path to hg38 not found

**Solution:** Update `DEFAULT_HG38_PATH` constant or omit `-hg` flag

---

## Testing

Test suite available in `test/phase1_testing/`

Run tests:
```bash
# Basic primer design
python phagpcr.py -f data/targets/CDS_0033_primase.fna \
  -o test/test1 -r 1 -nb 3

# With BLAST screening
python phagpcr.py -f data/targets/CDS_0033_primase.fna \
  -o test/test2 -r 1 -nb 3

# Gene extractor
python bz_gene_extractor.py -i data/genbank/bc14_KZ.gbk \
  -w "primase" -o test/test3 -sk 500
```

See `test/phase1_testing/TEST_RESULTS.md` for full test results.

---

## Project Status

| Component | Status | Notes |
|-----------|--------|-------|
| Workflow 1 (Target qPCR) | ✅ Complete | Fully tested and functional |
| Workflow 2 (Cocktail) | ✅ Complete | Includes dimer analysis |
| Workflow 3 (Unique regions) | 📋 Planned | See Phase 3 planning doc |
| Gene Extractor | ✅ Complete | Fully functional |
| Documentation | ✅ Complete | Updated 2025-11-13 |
| Tests | ✅ Passing | 4/4 test cases pass |

---

## Citation

If you use PHAGPCR in your research, please cite:

```
[Citation information to be added]
```

---

## License

GPL-3.0 License - see [LICENSE](LICENSE) file for details

---

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Follow PEP 8 coding standards (80 char line limit)
4. Add tests for new features
5. Submit a pull request

---

## Contact

For questions, issues, or feature requests, please open an issue on GitHub.

---

**Last Updated:** 2025-11-13
**Version:** 1.0.0
