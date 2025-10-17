# CHARLA Testing Suite

This directory contains comprehensive tests for the CHARLA Nextflow pipeline using nf-test.

## Directory Structure

```
tests/
├── README.md                    # This file
├── data/                        # Test data fixtures
│   ├── test_sample.fastq       # Small FASTQ for testing
│   ├── test_sample.fasta       # Small FASTA for testing
│   └── kmc_dbs/                # Mini KMC databases for testing (TODO)
├── modules/                     # Module-level tests
│   └── local/
│       ├── fastq_to_fasta.nf.test
│       ├── simplify_headers.nf.test
│       ├── generate_readmers_kmc.nf.test
│       ├── get_counts_kmc.nf.test
│       ├── generate_histogram_kmc.nf.test
│       └── ... (more module tests)
├── workflows/                   # Workflow-level tests (TODO)
│   └── charla_nf.nf.test
├── pipelines/                   # Full pipeline integration tests
│   └── default.nf.test         # End-to-end test
└── nextflow.config              # Test-specific config

```

## Running Tests

### Run all tests
```bash
nf-test test
```

### Run specific module tests
```bash
nf-test test tests/modules/local/fastq_to_fasta.nf.test
```

### Run tests with specific profile
```bash
nf-test test --profile docker
```

### Run tests in verbose mode
```bash
nf-test test --verbose
```

### Update snapshots
```bash
nf-test test --update-snapshot
```

## Test Categories

### 1. Module Tests (`modules/local/`)
Unit tests for individual processes:
- **fastq_to_fasta**: FASTQ → FASTA conversion with/without masking
- **simplify_headers**: FASTA header simplification
- **generate_readmers_kmc**: K-mer database generation
- **get_counts_kmc**: Per-read k-mer counting
- **generate_histogram_kmc**: K-mer frequency histograms
- **get_counts_cenhapmer_kmc**: Cenhapmer profiling (TODO)
- **combine_hybrid_profiles**: Hybrid read identification (TODO)
- **curate_hybrid_profiles**: Hybrid profile curation (TODO)
- **segment_reads**: Read segmentation by haplotype (TODO)
- **mapping_analysis**: Minimap2-based validation (TODO)

### 2. Workflow Tests (`workflows/`)
Integration tests for subworkflows:
- Complete preprocessing workflow
- K-mer profiling workflow
- Hybrid detection workflow

### 3. Pipeline Tests (`pipelines/`)
End-to-end tests:
- Full pipeline with minimal data
- Different tissue types (pollen/sperm)
- Different hybrid crosses

## Test Data

Test data is stored in `tests/data/` and includes:

### Minimal Test Files
- **test_sample.fastq**: 5 reads, ~250bp each, various quality profiles
- **test_sample.fasta**: 5 sequences with complex headers

### Full Test Dataset (requires download)
For integration testing, larger datasets are available:
```bash
# Download test datasets (TODO: add actual download location)
# wget -P tests/data/ https://...
```

## Writing New Tests

### Template for Module Test

```groovy
nextflow_process {
    name "Test Process YOUR_PROCESS"
    script "../your_process.nf"
    process "YOUR_PROCESS"
    tag "modules"
    tag "modules_local"
    tag "your_process"

    test("Should perform expected function") {
        when {
            process {
                """
                input[0] = [
                    [ id: 'test_id' ],
                    file('path/to/test/file')
                ]
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert process.out.your_output.size() == 1 },
                { assert snapshot(process.out).match() }
            )
        }
    }
}
```

### Best Practices

1. **Use descriptive test names**: Clearly state what the test validates
2. **Test edge cases**: Empty inputs, malformed data, boundary conditions
3. **Check outputs**: Verify file existence, format, and content where possible
4. **Use snapshots**: For complex outputs, use snapshots for regression testing
5. **Test failure modes**: Ensure processes fail gracefully with meaningful errors
6. **Use setup blocks**: For tests requiring upstream process outputs
7. **Tag appropriately**: Use tags for selective test running

## Continuous Integration

Tests are automatically run on:
- Pull requests
- Pushes to main/master
- Scheduled nightly builds

See `.github/workflows/ci.yml` for CI configuration.

## Troubleshooting

### Test Data Not Found
Ensure test data paths use `params.pipelines_testdata_base_path`:
```groovy
file(params.pipelines_testdata_base_path + 'charla_nf/test_sample.fastq', checkIfExists: true)
```

### KMC Tests Fail
KMC requires sufficient memory. Adjust in `nextflow.config`:
```groovy
process.memory = '8 GB'
```

### Snapshot Mismatches
After intentional changes, update snapshots:
```bash
nf-test test --update-snapshot
```

## Coverage

Current test coverage:
- ✅ FASTQ_TO_FASTA (5 tests)
- ✅ SIMPLIFY_HEADERS (5 tests)
- ✅ generate_readmers_kmc (5 tests)
- ✅ get_counts_kmc (2 tests)
- ✅ generate_histogram_kmc (2 tests)
- ⏳ get_counts_cenhapmer_kmc (pending)
- ⏳ combine_hybrid_profiles (pending)
- ⏳ curate_hybrid_profiles (pending)
- ⏳ segment_reads (pending)
- ⏳ mapping_analysis (pending)

## Contributing

When adding new modules or modifying existing ones:
1. Write tests first (TDD approach)
2. Ensure all tests pass before submitting PR
3. Aim for 100% module coverage
4. Document any special test requirements

## Resources

- [nf-test documentation](https://code.askimed.com/nf-test/)
- [Nextflow testing patterns](https://www.nextflow.io/docs/latest/testing.html)
- [nf-core testing guidelines](https://nf-co.re/docs/contributing/modules#writing-tests)
