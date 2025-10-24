# Rust Simulator - Technical Documentation

## Overview

The Rust simulator is the core sequence generation engine for the spacer benchmarking pipeline. It generates synthetic CRISPR spacer and contig sequences with controlled mutations, creating ground truth data for tool validation.

## Key Terminology

### CRITICAL: Understanding "Insertion"

The term "insertion" has **TWO DIFFERENT MEANINGS** in this codebase:

1. **Spacer Insertion (into contigs)**
   - **What it means**: Placing a spacer sequence INTO a contig at a random position
   - **Parameter**: `insertion_range` or `-ir` in CLI
   - **Example**: `-ir 1 3` means each spacer will be inserted 1-3 times into different contigs
   - **Purpose**: Creates the ground truth annotations (what we're trying to find)

2. **Insertion Mutation (indel)**
   - **What it means**: Adding extra bases WITHIN a spacer sequence as a mutation
   - **Parameter**: `n_insertion_range` or `-nir` in CLI  
   - **Example**: `-nir 0 2` means 0-2 insertion events within the spacer
   - **Purpose**: Simulates biological sequence variation (makes alignment harder)

Similarly for **Deletion Mutation**:
- **What it means**: Removing bases WITHIN a spacer sequence as a mutation
- **Parameter**: `n_deletion_range` or `-ndr` in CLI
- **Example**: `-ndr 0 1` means 0-1 deletion events within the spacer

## Simulation Parameters

### Required Parameters

| Parameter (Rust) | CLI Flag | Description | Example |
|-----------------|----------|-------------|---------|
| `contig_length_range` | `-cl` | Min and max contig lengths | `1000 5000` |
| `spacer_length_range` | `-ls` | Min and max spacer lengths | `20 50` |
| `sample_size_contigs` | `-nc` | Number of contigs to generate | `100` |
| `sample_size_spacers` | `-ns` | Number of spacers to generate | `50` |
| `insertion_range` | `-ir` | **How many times to insert each spacer** | `1 1` |
| `n_mismatch_range` | `-lm` | Substitution mutations per spacer | `0 0` |
| `n_insertion_range` | `-nir` | **Insertion mutations per spacer** | `0 0` |
| `n_deletion_range` | `-ndr` | **Deletion mutations per spacer** | `0 0` |
| `prop_rc` | `-prc` | Proportion of reverse complement | `0.5` |

### Optional Parameters

| Parameter (Rust) | CLI Flag | Description | Default |
|-----------------|----------|-------------|---------|
| `id_prefix` | `-id` | Prefix for sequence IDs | Auto-generated hash |
| `gc_content` | `--gc-content` | GC content percentage (0-100) | Random |
| `contig_distribution` | `--contig-distribution` | Length distribution type | `uniform` |
| `spacer_distribution` | `--spacer-distribution` | Length distribution type | `uniform` |
| `a_frac`, `t_frac`, etc. | `--a-frac`, etc. | Individual base fractions | Equal |

## Simulation Algorithm

### High-Level Workflow

```
1. Generate random contigs
   ↓
2. Generate random spacers
   ↓
3. For each spacer:
   a. Determine number of insertions (from insertion_range)
   b. For each insertion:
      - Select random contig
      - Select random position in contig
      - Determine if reverse complement (from prop_rc)
      - Determine mutations:
        * n_mismatches (substitutions)
        * n_insertions (insertion indels)
        * n_deletions (deletion indels)
      - Apply mutations to create spacer variant
      - Insert variant into contig
      - Record ground truth (position, strand, mutations)
   ↓
4. Output:
   - simulated_contigs.fa
   - simulated_spacers.fa
   - planned_ground_truth.tsv
```

### Detailed Implementation

#### 1. Contig Generation

```rust
// Generate random contig sequences
fn generate_random_sequence(
    length: usize,
    base_composition: &BaseComposition
) -> String
```

- Uses weighted random sampling based on base composition (GC content or individual fractions)
- Applies length distribution (uniform, normal, or bell curve)

#### 2. Spacer Generation

```rust
// Generate spacer sequences
// Similar to contig generation but with spacer length range
```

#### 3. Spacer Insertion with Mutations

**Key Code Section** (from `src/rust_simulator/src/lib.rs`):

```rust
// Determine number of insertions for this spacer
let n_insertions = rng.random_range(insertion_range.0..=insertion_range.1);

// Create insertion plans
let mut plans = Vec::with_capacity(n_insertions);
for _ in 0..n_insertions {
    let is_rc = rng.random_bool(prop_rc);
    let n_mismatches = rng.random_range(n_mismatch_range.0..=n_mismatch_range.1);
    let n_insertions = rng.random_range(n_insertion_range.0..=n_insertion_range.1);
    let n_deletions = rng.random_range(n_deletion_range.0..=n_deletion_range.1);
    plans.push((is_rc, n_mismatches, n_insertions, n_deletions));
}
```

**Understanding this code:**
- `n_insertions` (first line): How many TIMES to insert the spacer
- `n_insertions` (in loop): How many insertion MUTATIONS to apply
- Each insertion plan includes: strand orientation, substitutions, insertion indels, deletion indels

#### 4. Mutation Application

```rust
fn apply_mutations(
    sequence: &str,
    n_mismatches: usize,
    n_insertions: usize,  // insertion indels
    n_deletions: usize    // deletion indels
) -> String
```

**Mutation Order:**
1. **Substitutions**: Replace random bases with different bases
2. **Insertions**: Add random bases at random positions
3. **Deletions**: Remove bases at random positions

**Important Notes:**
- Mutations are applied to the spacer BEFORE insertion into contig
- The original spacer sequence is preserved in `simulated_spacers.fa`
- The mutated variant is what gets inserted into the contig
- Ground truth records the ACTUAL inserted position and length

#### 5. Position Calculation

```rust
// Insert spacer variant into contig
let start_pos = /* random position */;
let actual_end_pos = start_pos + spacer_variant.len();
```

**Critical Detail:**
- Start position is 0-based (Python-style indexing)
- End position is 0-based EXCLUSIVE (Python slice style)
- Length can differ from original spacer if indels were applied
- Strand is recorded as boolean: `true` = reverse complement, `false` = forward

## Ground Truth Format

### Output File: `planned_ground_truth.tsv`

| Column | Type | Description |
|--------|------|-------------|
| `spacer_id` | String | Original spacer identifier |
| `contig_id` | String | Contig where spacer was inserted |
| `start` | Integer | Start position (0-based, inclusive) |
| `end` | Integer | End position (0-based, exclusive) |
| `strand` | Boolean | `true` = reverse complement, `false` = forward |
| `mismatches` | Integer | Number of substitution mutations applied |

**Note**: The `mismatches` field only counts substitutions, NOT indels.

## Common Use Cases

### 1. Perfect Match Validation

Test alignment tools with zero mutations:

```bash
spacer_bencher simulate \
  -nc 100 -ns 50 \
  -ir 1 1 \    # Insert each spacer exactly once
  -lm 0 0 \    # No substitutions
  -nir 0 0 \   # No insertion indels
  -ndr 0 0 \   # No deletion indels
  -o tests/perfect
```

**Expected result**: Tools should find 100% exact matches.

### 2. Substitution-Only Testing

Test tools with only substitution mutations:

```bash
spacer_bencher simulate \
  -nc 100 -ns 50 \
  -ir 1 1 \    # Single insertion per spacer
  -lm 0 3 \    # 0-3 substitutions
  -nir 0 0 \   # No insertion indels
  -ndr 0 0 \   # No deletion indels
  -o tests/mismatches_only
```

**Expected result**: All matches should have same length as original spacer.

### 3. Full Mutation Spectrum

Test tools with all mutation types:

```bash
spacer_bencher simulate \
  -nc 500 -ns 100 \
  -ir 1 3 \    # 1-3 insertions per spacer
  -lm 0 2 \    # 0-2 substitutions
  -nir 0 1 \   # 0-1 insertion indels
  -ndr 0 1 \   # 0-1 deletion indels
  -prc 0.5 \   # 50% reverse complement
  -o tests/full_mutations
```

**Expected result**: Variable-length alignments with mixed mutation types.

## Performance Considerations

### Speed Optimizations

1. **Parallel Processing**: The simulator uses Rayon for parallel contig/spacer generation
   ```rust
   use rayon::prelude::*;
   contigs.par_iter_mut().for_each(|contig| { ... });
   ```

2. **Pre-allocation**: Vectors are pre-allocated to avoid repeated allocations
   ```rust
   let mut plans = Vec::with_capacity(n_insertions);
   ```

3. **Efficient String Operations**: Uses byte-level operations where possible

### Memory Usage

- Typical memory usage: ~50-100 MB for standard simulations
- Scales linearly with:
  - Number of contigs × average contig length
  - Number of spacers × average spacer length × insertions per spacer

### Typical Runtime

On a modern workstation (8 threads):
- Small (100 contigs, 50 spacers): ~1-2 seconds
- Medium (1000 contigs, 200 spacers): ~5-10 seconds
- Large (10000 contigs, 1000 spacers): ~30-60 seconds

## Validation and Debugging

### Verify Mode

The simulator includes a verification mode that checks:

```rust
// Verify that inserted spacers can be found in contigs
for ground_truth_entry in ground_truth {
    let contig = &contigs[ground_truth_entry.contig_id];
    let region = &contig[ground_truth_entry.start..ground_truth_entry.end];
    // Verify region matches expected spacer variant
}
```

Enable with `--verify` flag:

```bash
spacer_bencher simulate --verify -o tests/validated
```

### Common Issues

1. **Empty ground truth**: 
   - Check that `insertion_range` is not `(0, 0)`
   - Ensure contigs are long enough to accommodate spacers

2. **Unexpected indels**:
   - Remember: `-nir` and `-ndr` create indels, `-ir` controls placement count

3. **Length mismatches**:
   - Ground truth `end - start` includes indels
   - Original spacer length is in `simulated_spacers.fa`

## Python Bindings

The Rust simulator is exposed to Python via PyO3:

```python
from rust_simulator import simulate_data

contigs, spacers, ground_truth = simulate_data(
    contig_length_range=(1000, 5000),
    spacer_length_range=(20, 50),
    # ... other parameters
)
```

## Future Enhancements

Potential improvements:

1. **Non-uniform mutation rates**: Allow specifying different mutation rates per spacer
2. **Clustered mutations**: Simulate mutation hotspots
3. **Structured variations**: Support specific mutation patterns (e.g., always at ends)
4. **Quality scores**: Generate mock quality scores for more realistic data

## Contributing

When modifying the Rust simulator:

1. Update this documentation
2. Add tests for new features
3. Benchmark performance changes
4. Update Python bindings if interface changes

## References

- Rust implementation: `src/rust_simulator/src/lib.rs`
- Python wrapper: `src/bench/utils/functions.py` (`simulate_data_rust()`)
- CLI interface: `src/bench/cli.py` (`simulate` command)
