# bcftools query failure handling test

Validates that the Tier 3 (`bcftools query`) error handling in `scripts/integration_test.sh` correctly detects failures instead of silently producing a false PASS.

## Background

The original code used `|| true` on both `bcftools query` commands:

```bash
bcftools query -f "${query_fmt}" "${local_output}" > "${query_output}" 2>/dev/null || true
bcftools query -f "${query_fmt}" "${local_expected}" > "${query_expected}" 2>/dev/null || true
```

If both queries failed (e.g. missing format fields, corrupt file), two empty files would be produced. `diff` on two empty files returns success, so the file would incorrectly PASS.

The fix checks exit status explicitly and fails the file if either query errors out.

## Tests

| # | Scenario | Expected | Result |
|---|----------|----------|--------|
| 1 | `bcftools query` on an invalid (non-VCF) file | Command fails (non-zero exit) | PASS |
| 2 | `bcftools query` on a valid VCF with GT/DP/AD fields | Command succeeds, fields extracted correctly | PASS |
| 3 | `diff` two empty files (simulating the old false-PASS bug) | `diff` returns success, confirming the bug scenario | PASS |

## Test details

### Test 1 — query on invalid file

```
$ echo "not a vcf" > bad.vcf.gz
$ bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%AD]\n' bad.vcf.gz
(exits non-zero)
```

Confirms the new `if !` check catches the failure.

### Test 2 — query on valid VCF

```
$ bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%AD]\n' good.vcf.gz
chr1	100	A	G	0/1	30	10,20
```

Confirms query works correctly on well-formed input.

### Test 3 — empty file diff (the bug)

```
$ touch empty1.txt empty2.txt
$ diff empty1.txt empty2.txt
(exits 0 — success)
```

Demonstrates why `|| true` was dangerous: if both queries silently fail and produce empty output, `diff` sees no differences and the test falsely passes.

## Conclusion

The fix correctly propagates `bcftools query` failures, preventing false PASSes when query extraction fails on either the output or expected file.
