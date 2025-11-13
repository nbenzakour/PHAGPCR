# Workflow 3 (UNIQUE TARGET) Implementation Plan

**Status:** Planned - Not Yet Implemented
**Target:** runtype=3
**Date:** 2025-11-13

---

## Objective

Implement unique region identification for phage genomes when specific genes (like primase) don't have unique primer-binding sites. Find unique regions anywhere in the phage genome that can distinguish it from a background database.

---

## Use Case

- **Input:** Phage genome (draft or complete) + Background database (other phages/bacteria)
- **Goal:** Design primers that will ONLY amplify the target phage, not the background
- **Application:** Diagnostic/detection when conserved genes cannot provide specificity

---

## Proposed Algorithm

```
1. Input: Target phage genome + Background database (other phages/bacteria)
2. Identify unique regions in target genome vs database
3. Filter regions by size (must accommodate primers + amplicon)
4. Run Primer3 on unique regions
5. Screen primers with MFEprimer (bacteria DB + human genome)
6. Output: Specific primers for unpredictable/unique target regions
```

---

## Implementation Complexity Breakdown

### Component 1: Unique Region Finder ⚠️ MODERATE-HIGH

**Approach: BLAST-based (Recommended)**

```python
# Pseudocode
1. Split target genome into sliding windows (e.g., 500bp, step 250bp)
2. BLAST each window against background database
3. Windows with no hits (or only self-hits) = unique regions
4. Merge overlapping unique windows
5. Filter by minimum length (e.g., >200bp for primer design)
```

**Implementation Details:**
- BLAST subprocess calls (similar to existing MFEprimer code)
- Window generation logic: ~50 lines
- BLAST result parsing: ~100 lines
- Region merging: ~50 lines

**Estimated Complexity:** 200-300 lines of code

**Alternative Approach: K-mer based (Faster but less sensitive)**
```python
# Pseudocode
1. Build k-mer index of background database (e.g., k=31)
2. Slide k-mers across target genome
3. Regions with unique k-mers = candidate unique regions
4. Validate with BLAST for confirmation
```

**Estimated Complexity:** 250-350 lines

---

### Component 2: Integration with Primer3 ✅ LOW

**What's needed:**
- Extract unique region sequences
- Pass to existing `get_primers()` function
- Already mostly implemented!

**Estimated Complexity:** 20-30 lines

---

### Component 3: Screening Pipeline ✅ LOW

**What's needed:**
- Use existing MFEprimer functions
- Screen against bacteria DB only (not phage DB - would defeat purpose)
- Screen against human genome

**Estimated Complexity:** Already implemented, just needs calling (~10 lines)

---

### Component 4: Output & Reporting ⚠️ MODERATE

**What's needed:**
- Report which regions are unique
- Show genomic coordinates
- Indicate confidence (how unique is unique?)
- Summary statistics

**Estimated Complexity:** 50-100 lines

---

## Total Implementation Estimate

| Component | Lines of Code | Complexity | Time Estimate |
|-----------|--------------|------------|---------------|
| Unique region finder (BLAST) | 250-300 | Moderate-High | 3-4 hours |
| Primer3 integration | 20-30 | Low | 30 min |
| MFEprimer screening | 10 | Low | 15 min |
| Output & reporting | 50-100 | Moderate | 1-2 hours |
| Testing & debugging | - | Moderate | 2-3 hours |
| **TOTAL** | **~330-440 lines** | **Moderate** | **7-10 hours** |

---

## Implementation Plan

### Phase A: Core Functionality (MVP)

#### 1. BLAST-based unique region finder
- Sliding window approach
- BLAST against background DB
- Identify regions with no hits (or below threshold)

#### 2. Extract and design primers
- Use existing `get_primers()` function
- Filter by Primer3 success

#### 3. Basic screening
- MFEprimer against bacteria-only DB
- MFEprimer against human genome

#### 4. Simple output
- CSV with unique regions (coordinates, length, uniqueness score)
- CSV with designed primers
- Summary statistics

**Estimated effort:** 6-8 hours of focused work

---

### Phase B: Enhancements (Optional)

#### 5. Advanced filtering
- GC content filtering (avoid extreme GC%)
- Repeat region filtering
- Low complexity filtering

#### 6. Confidence scoring
- How unique is the region (0 hits vs 1-2 hits vs many)
- Phylogenetic distance to nearest hit
- Uniqueness score metric

#### 7. Visualization
- Genome map showing unique regions
- Primer locations on genome
- Coverage plots

**Additional effort:** 4-6 hours

---

## Key Design Decisions

### 1. Window size for uniqueness detection?
- **Recommendation:** 500bp windows with 250bp overlap
- **Rationale:** Must be large enough for amplicon (100-300bp) + primers (40bp)
- **Configurable:** Yes, add as parameter

### 2. What defines "unique"?
**Options:**
- Zero BLAST hits? (very strict, may find nothing)
- No hits above certain identity threshold (e.g., <80%)? (more permissive)
- No hits above coverage threshold (e.g., <50% query coverage)?

**Recommendation:** Configurable, default <80% identity over <50% query coverage

### 3. Minimum unique region size?
- **Recommendation:** 200bp minimum
- **Rationale:** Allows for primer design + small amplicon
- **Configurable:** Yes, add as parameter `--min-unique-length`

### 4. How many unique regions to report?
**Options:**
- All unique regions?
- Top N by uniqueness score?
- Best scoring only?

**Recommendation:** Top 10 by uniqueness score, but save all in detailed output

---

## Command Line Interface

### Proposed Arguments for runtype=3

```bash
python phagpcr.py \
  -f phage_genome.fasta \           # Target phage genome
  -o results/ \                      # Output directory
  -b bacteria_background.fna \       # Background database to check against
  -r 3 \                             # Runtype 3: unique region detection
  -hg \                              # Also screen against human genome
  --window-size 500 \                # Sliding window size (bp)
  --window-step 250 \                # Window step size (bp)
  --min-unique-length 200 \          # Minimum unique region length (bp)
  --identity-threshold 80 \          # Max % identity to background (below = unique)
  --coverage-threshold 50 \          # Max % coverage to background (below = unique)
  -nb 3                              # Number of primer pairs per unique region
```

---

## Implementation Steps

### Step 1: Create unique region finder function
```python
def find_unique_regions(genome_file, background_db, window_size=500,
                       window_step=250, identity_cutoff=80,
                       coverage_cutoff=50, min_length=200):
    """
    Identify unique regions in genome vs background database.

    Args:
        genome_file: Path to target genome FASTA
        background_db: Path to background BLAST database
        window_size: Size of sliding window (bp)
        window_step: Step size for sliding window (bp)
        identity_cutoff: Max % identity to consider unique
        coverage_cutoff: Max % coverage to consider unique
        min_length: Minimum unique region length (bp)

    Returns:
        DataFrame with columns: region_id, start, end, length,
                                uniqueness_score, sequence
    """
    # Implementation here
```

### Step 2: BLAST wrapper function
```python
def blast_window_against_db(window_seq, db_path, temp_dir):
    """
    BLAST a sequence window against database.

    Returns:
        dict with best hit info (identity, coverage, eval) or None if no hits
    """
    # Implementation here
```

### Step 3: Region merging function
```python
def merge_overlapping_regions(regions_df, min_length):
    """
    Merge overlapping unique windows into contiguous regions.

    Returns:
        DataFrame with merged regions
    """
    # Implementation here
```

### Step 4: Integration with main workflow
```python
elif runtype == 3:
    print("Step 1 - Identify unique regions in genome")
    unique_regions = find_unique_regions(fasta_file, blast_db, ...)

    print("Step 2 - Design primers for unique regions")
    for region in unique_regions:
        primers = get_primers(region.id, region.sequence, primer3_params)
        # ... rest of workflow

    print("Step 3 - Screen primers with MFEprimer")
    # Use existing screening functions
```

---

## Testing Strategy

### Test Case 1: Known unique region
- Use phage with characterized unique region
- Verify algorithm identifies it correctly

### Test Case 2: Highly conserved phage
- Use phage very similar to database
- Verify it finds the few unique regions (or reports none)

### Test Case 3: Novel phage
- Use phage very different from database
- Verify it finds many unique regions

### Test Case 4: Edge cases
- Very small genome (<5kb)
- Very large genome (>200kb)
- Circular vs linear genomes

---

## Output Files

### For runtype=3, generate:

1. **unique_regions.csv**
   - Columns: region_id, start, end, length, gc_content, uniqueness_score,
             num_hits, best_hit_identity, best_hit_coverage

2. **unique_regions.fna**
   - FASTA file with all unique region sequences

3. **primers_from_unique_regions.csv**
   - Same format as existing prefiltered_primers_file.csv
   - Added columns: source_region_id, region_start, region_end

4. **mfe_screening_results.csv**
   - MFEprimer screening against bacteria DB

5. **mfe_spec_results_hg38.csv** (if -hg flag)
   - MFEprimer screening against human genome

6. **summary_report.txt**
   - Number of unique regions found
   - Total unique sequence length
   - Number of primer pairs designed
   - Number passing screening

---

## Dependencies

### New Dependencies Required:
- ✅ BLAST+ (blastn) - likely already installed
- ✅ BioPython - already in requirements.txt
- ✅ pandas - already in requirements.txt

### No new Python packages needed!

---

## Potential Issues & Solutions

### Issue 1: BLAST is slow for large databases
**Solution:**
- Use BLAST with `-max_target_seqs 10` to speed up
- Consider pre-filtering with k-mers (Phase B enhancement)
- Cache BLAST results

### Issue 2: No unique regions found
**Solution:**
- Relax identity/coverage thresholds
- Reduce minimum region length
- Report best "least conserved" regions even if not truly unique

### Issue 3: Too many unique regions
**Solution:**
- Rank by uniqueness score
- Prioritize longer regions
- Prioritize regions with better GC content for primer design

### Issue 4: Circular genome handling
**Solution:**
- Detect circularity from FASTA header or sequence features
- Handle wrap-around windows specially
- Report coordinates relative to start position

---

## Success Criteria

Workflow 3 implementation is complete when:

- ✅ Can identify unique regions in target genome vs background database
- ✅ Can design primers on unique regions using Primer3
- ✅ Can screen primers with MFEprimer
- ✅ Generates all required output files
- ✅ Passes test cases
- ✅ Documentation updated (README, help text)
- ✅ Example workflow in documentation

---

## Timeline

### Recommended Phased Approach:

**Phase 1 & 2 (Code Quality):** ~3-4 hours
- Complete type hints, logging, pathlib, README updates
- Merge to main

**Phase 3A (Core Unique Region Detection):** ~6-8 hours
- Implement BLAST-based unique region finder
- Integrate with Primer3
- Basic MFEprimer screening
- Test with sample data

**Phase 3B (Enhancements - Optional):** ~4-6 hours
- Advanced filtering
- Confidence scoring
- Visualization

**Total for complete Phase 1+2+3A:** ~10-12 hours
**Total for complete Phase 1+2+3A+3B:** ~15-20 hours

---

## Open Questions

1. Should we support multiple target genomes simultaneously?
2. How to handle draft genomes with multiple contigs?
3. Should we integrate with existing phage genome databases (INPHARED, etc.)?
4. Do we want to parallelize BLAST searches for speed?
5. Should unique regions be validated experimentally before deployment?

---

**Last Updated:** 2025-11-13
**Status:** Planning complete, ready for implementation
