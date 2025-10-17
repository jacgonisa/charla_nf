# CHARLA Pipeline Generalization Plan

## Current Hardcoded Assumptions

### 1. **Accession Names (Col/Ler)**

#### Locations in Code:
- **`main.nf:60-67`**: Hardcoded `['Col', 'Ler']` for k-mer combinations
- **`workflows/charla_nf.nf:392-405`**: Parameters `--col_bed_file` and `--ler_bed_file`
- **`workflows/charla_nf.nf:410-422`**: PAF file filtering uses "Col" and "Ler" in filenames
- **`nextflow.config:14`**: Default example uses `ColLer_F1_pollen.fq.gz`
- **`nextflow.config:43-46`**: Parameters named `col_bed_file` and `ler_bed_file`
- **`README.md`**: All examples use Col×Ler cross

#### Scripts that may contain hardcoded names:
- Rust binaries may have hardcoded accession assumptions
- Plotting scripts likely expect "Col" and "Ler" in output names

### 2. **Chromosome Number (1-5)**

#### Locations:
- **`main.nf:60-68`**: Range `(1..5)` hardcoded for Arabidopsis chromosomes
- May affect scripts that parse or filter by chromosome

### 3. **File Naming Conventions**

#### Patterns that assume Col/Ler:
- `againstCol.paf` / `againstLer.paf`
- BED files named with Col/Ler prefixes
- Reference genome directory structure may assume Col/Ler naming

## Proposed Generalization Strategy

### Phase 1: Parameter Abstraction

Replace accession-specific parameters with generic parent names:

```groovy
// OLD (Col/Ler specific):
params {
    col_bed_file = null
    ler_bed_file = null
}

// NEW (generic):
params {
    parent1_name = 'Parent1'  // e.g., 'Col'
    parent2_name = 'Parent2'  // e.g., 'Ler'
    parent1_bed = null
    parent2_bed = null
    num_chromosomes = 5  // Default for Arabidopsis
}
```

### Phase 2: Dynamic Channel Creation

Update `main.nf` to create k-mer combinations dynamically:

```groovy
// OLD:
ch_kmc_combinations = Channel
    .from(['Col', 'Ler'])
    .flatMap { ... }

// NEW:
ch_kmc_combinations = Channel
    .from([params.parent1_name, params.parent2_name])
    .flatMap { acc ->
        ['CEN', 'ARMS'].collectMany { cat ->
            (1..params.num_chromosomes).collect { chr ->
                [acc, cat, chr]
            }
        }
    }
```

### Phase 3: Sample Sheet Schema

Create a sample sheet for batch processing multiple samples:

```csv
sample_id,reads,input_format,tissue_type,cross,parent1,parent2
ColLer_F1_pollen,data/pollen.fq.gz,fastq,pollen,ColxLer,Col,Ler
ColLer_F1_sperm,data/sperm.fq.gz,fastq,sperm,ColxLer,Col,Ler
```

Schema definition (`assets/schema_input.json`):

```json
{
    "$schema": "http://json-schema.org/draft-07/schema",
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "sample_id": {"type": "string", "pattern": "^\\S+$"},
            "reads": {"type": "string", "format": "file-path"},
            "input_format": {"type": "string", "enum": ["fastq", "fasta"]},
            "tissue_type": {"type": "string", "enum": ["pollen", "sperm"]},
            "cross": {"type": "string"},
            "parent1": {"type": "string"},
            "parent2": {"type": "string"}
        },
        "required": ["sample_id", "reads", "input_format", "parent1", "parent2"]
    }
}
```

### Phase 4: Update Module Inputs

Modules that need parent names:
- `MAPPING_ANALYSIS` - needs to know parent names for PAF filtering
- `PLOT_CROSSOVER_MAP` - needs parent BED files and names
- `NON_HYBRID_READS_ANALYSIS` - may need parent info

### Phase 5: Backward Compatibility

Keep old parameter names as aliases (deprecated):

```groovy
// Backward compatibility
if (params.col_bed_file && !params.parent1_bed) {
    params.parent1_bed = params.col_bed_file
    params.parent1_name = 'Col'
    log.warn "⚠️  DEPRECATED: --col_bed_file is deprecated. Use --parent1_bed instead."
}
```

### Phase 6: Tissue-Type Specific Logic

Add optional tissue-type metadata for future differential processing:

```groovy
params {
    tissue_type = null  // Optional: 'pollen' or 'sperm'
}
```

This allows for tissue-specific filtering or analysis in future versions.

## Implementation Order

1. ✅ Create this planning document
2. ⏳ Update `nextflow.config` with new parameters
3. ⏳ Update `main.nf` k-mer combination logic
4. ⏳ Update `workflows/charla_nf.nf` BED file handling
5. ⏳ Create sample sheet schema
6. ⏳ Add backward compatibility layer
7. ⏳ Update README with generalized examples
8. ⏳ Test with different hybrid cross data

## Testing Strategy

Test cases:
1. **Backward compatibility**: Run with old `--col_bed_file` / `--ler_bed_file` params
2. **New parameters**: Run with `--parent1_name` / `--parent2_name`
3. **Different species**: Test with 7-chromosome organism (if available)
4. **Sample sheet**: Test batch processing with multiple samples

## Benefits

1. **Flexibility**: Works with any F1 hybrid cross, not just Col×Ler
2. **Scalability**: Sample sheet enables batch processing
3. **Clarity**: Generic parent1/parent2 naming is more intuitive
4. **Future-proof**: Easy to add tissue-specific processing later
5. **Reusability**: Can be applied to other species with different chromosome numbers
