#!/usr/bin/env bash
#
# Integration test: invoke the VCF normalisation Lambda on test VCFs and
# compare outputs against expected files using three levels of comparison:
#   1. bcftools stats  — summary statistics sanity check
#   2. bcftools isec   — site-level (CHROM/POS/REF/ALT) comparison
#   3. bcftools query  — field-level (GT/DP/AD) comparison
#
# A markdown report is written to $REPORT_FILE (default: integration_report.md).
#
# Usage:
#   ./scripts/integration_test.sh <bucket> [input_prefix] [expected_prefix]
#
# Environment variable overrides:
#   FUNCTION_NAME  Lambda function name        (default: vcf-normalisation)
#   MAX_PARALLEL   Concurrent Lambda invokes   (default: 10)
#   POLL_TIMEOUT   Per-file poll timeout in s   (default: 120)
#   REPORT_FILE    Markdown report path         (default: integration_report.md)

set -euo pipefail

# ── Arguments & defaults ─────────────────────────────────────────────────────

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <bucket> [input_prefix] [expected_prefix]" >&2
    exit 1
fi

BUCKET="$1"
INPUT_PREFIX="${2:-test/input/}"
EXPECTED_PREFIX="${3:-test/expected/}"
FUNCTION_NAME="${FUNCTION_NAME:-vcf-normalisation}"
MAX_PARALLEL="${MAX_PARALLEL:-10}"
POLL_TIMEOUT="${POLL_TIMEOUT:-120}"
REPORT_FILE="${REPORT_FILE:-integration_report.md}"

# Derive output prefix by replacing the last /input/ with /output/ (matches
# the Lambda's rfind-based replacement logic)
OUTPUT_PREFIX=$(echo "${INPUT_PREFIX}" | sed 's|\(.*\)/input/|\1/output/|')

# ── Discover test VCFs ───────────────────────────────────────────────────────

mapfile -t FILES < <(
    aws s3 ls "s3://${BUCKET}/${INPUT_PREFIX}" \
        | awk '{print $NF}' \
        | grep '\.vcf\.gz$'
)

if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "No .vcf.gz files found in s3://${BUCKET}/${INPUT_PREFIX}" >&2
    exit 1
fi

echo "Found ${#FILES[@]} test VCFs in s3://${BUCKET}/${INPUT_PREFIX}"

# ── Clean previous outputs ───────────────────────────────────────────────────

echo "Cleaning previous outputs in s3://${BUCKET}/${OUTPUT_PREFIX}..."
aws s3 rm "s3://${BUCKET}/${OUTPUT_PREFIX}" --recursive > /dev/null 2>&1 || true

# ── Invoke Lambda for each file ─────────────────────────────────────────────

echo "Invoking Lambda for ${#FILES[@]} files (max ${MAX_PARALLEL} parallel)..."

invoke_lambda() {
    local file="$1"
    local key="${INPUT_PREFIX}${file}"
    local payload
    payload=$(jq -nc --arg b "${BUCKET}" --arg k "${key}" '{bucket: $b, key: $k}')

    aws lambda invoke \
        --function-name "${FUNCTION_NAME}" \
        --invocation-type Event \
        --payload "${payload}" \
        --cli-binary-format raw-in-base64-out \
        /dev/null > /dev/null 2>&1
}

sent=0
total=${#FILES[@]}
for file in "${FILES[@]}"; do
    invoke_lambda "${file}" &
    sent=$((sent + 1))

    # Throttle to MAX_PARALLEL concurrent jobs
    if (( sent % MAX_PARALLEL == 0 )); then
        wait
        printf "  [%d/%d] invocations sent\r" "${sent}" "${total}"
    fi
done
wait
echo "  [${total}/${total}] invocations sent"

# ── Wait for outputs to appear ───────────────────────────────────────────────

echo "Waiting for outputs in s3://${BUCKET}/${OUTPUT_PREFIX}..."

ready=0
missing=()
for file in "${FILES[@]}"; do
    output_key="${OUTPUT_PREFIX}${file}"
    elapsed=0
    found=true

    while ! aws s3 ls "s3://${BUCKET}/${output_key}" > /dev/null 2>&1; do
        if (( elapsed >= POLL_TIMEOUT )); then
            found=false
            missing+=("${file}")
            break
        fi
        sleep 5
        elapsed=$((elapsed + 5))
    done

    if [[ "${found}" == true ]]; then
        ready=$((ready + 1))
    fi
    printf "  [%d/%d] outputs ready\r" "${ready}" "${total}"
done
echo "  [${ready}/${total}] outputs ready"

if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: ${#missing[@]} files timed out waiting for output:" >&2
    for f in "${missing[@]}"; do
        echo "  - ${f}" >&2
    done
fi

# ── Compare results ──────────────────────────────────────────────────────────

echo "Comparing results..."

TMPDIR=$(mktemp -d)
trap 'rm -rf "${TMPDIR}"' EXIT

passed=0
failed=0
failures=()

# Per-file result arrays for the report
result_status=()
result_reason=()
result_input_count=()
result_expected_count=()
result_output_count=()
result_detail=()
result_full_diff=()
result_isec_output_only=()
result_isec_expected_only=()

for i in "${!FILES[@]}"; do
    file="${FILES[$i]}"

    # Skip files that never appeared
    if printf '%s\n' "${missing[@]}" 2>/dev/null | grep -qxF "${file}"; then
        echo "  [SKIP] ${file} — output timed out"
        failed=$((failed + 1))
        failures+=("${file}: output timed out")
        result_status[$i]="FAIL"
        result_reason[$i]="output timed out"
        result_input_count[$i]=""
        result_expected_count[$i]=""
        result_output_count[$i]=""
        result_detail[$i]=""
        result_full_diff[$i]=""
        result_isec_output_only[$i]=""
        result_isec_expected_only[$i]=""
        continue
    fi

    input_key="${INPUT_PREFIX}${file}"
    output_key="${OUTPUT_PREFIX}${file}"
    expected_key="${EXPECTED_PREFIX}${file}"

    local_input="${TMPDIR}/input_${file}"
    local_output="${TMPDIR}/output_${file}"
    local_expected="${TMPDIR}/expected_${file}"
    isec_dir="${TMPDIR}/isec_${file%.vcf.gz}"

    # Download input, output, and expected
    if ! aws s3 cp "s3://${BUCKET}/${input_key}" "${local_input}" > /dev/null 2>&1; then
        echo "  [FAIL] ${file} — could not download input"
        failed=$((failed + 1))
        failures+=("${file}: could not download input")
        result_status[$i]="FAIL"
        result_reason[$i]="could not download input"
        result_input_count[$i]=""
        result_expected_count[$i]=""
        result_output_count[$i]=""
        result_detail[$i]=""
        result_full_diff[$i]=""
        result_isec_output_only[$i]=""
        result_isec_expected_only[$i]=""
        continue
    fi

    if ! aws s3 cp "s3://${BUCKET}/${output_key}" "${local_output}" > /dev/null 2>&1; then
        echo "  [FAIL] ${file} — could not download output"
        failed=$((failed + 1))
        failures+=("${file}: could not download output")
        result_status[$i]="FAIL"
        result_reason[$i]="could not download output"
        result_input_count[$i]=""
        result_expected_count[$i]=""
        result_output_count[$i]=""
        result_detail[$i]=""
        result_full_diff[$i]=""
        result_isec_output_only[$i]=""
        result_isec_expected_only[$i]=""
        continue
    fi

    if ! aws s3 cp "s3://${BUCKET}/${expected_key}" "${local_expected}" > /dev/null 2>&1; then
        echo "  [FAIL] ${file} — could not download expected"
        failed=$((failed + 1))
        failures+=("${file}: could not download expected")
        result_status[$i]="FAIL"
        result_reason[$i]="could not download expected"
        result_input_count[$i]=""
        result_expected_count[$i]=""
        result_output_count[$i]=""
        result_detail[$i]=""
        result_full_diff[$i]=""
        result_isec_output_only[$i]=""
        result_isec_expected_only[$i]=""
        continue
    fi

    # Index all files
    bcftools index "${local_input}" 2>/dev/null || true
    bcftools index "${local_output}" 2>/dev/null || true
    bcftools index "${local_expected}" 2>/dev/null || true

    # Count variants
    input_count=$(bcftools view -H "${local_input}" 2>/dev/null | wc -l)
    output_count=$(bcftools view -H "${local_output}" 2>/dev/null | wc -l)
    expected_count=$(bcftools view -H "${local_expected}" 2>/dev/null | wc -l)
    result_input_count[$i]="${input_count}"
    result_expected_count[$i]="${expected_count}"
    result_output_count[$i]="${output_count}"

    # ── Tier 1: bcftools stats ──────────────────────────────────────────────

    stats_file="${TMPDIR}/stats_${file%.vcf.gz}.txt"
    if ! bcftools stats "${local_output}" "${local_expected}" > "${stats_file}" 2>/dev/null; then
        echo "  [FAIL] ${file} — bcftools stats failed"
        failed=$((failed + 1))
        failures+=("${file}: bcftools stats failed")
        result_status[$i]="FAIL"
        result_reason[$i]="bcftools stats failed"
        result_detail[$i]=""
        result_full_diff[$i]=""
        result_isec_output_only[$i]=""
        result_isec_expected_only[$i]=""
        continue
    fi

    # Parse SN lines — compare record counts between the two files
    output_records=$(awk -F'\t' '/^SN\t0\tnumber of records:/ {print $4}' "${stats_file}")
    expected_records=$(awk -F'\t' '/^SN\t1\tnumber of records:/ {print $4}' "${stats_file}")

    if [[ "${output_records}" != "${expected_records}" ]]; then
        echo "  [FAIL] ${file} — record count mismatch (output: ${output_records}, expected: ${expected_records})"
        failed=$((failed + 1))
        failures+=("${file}: record count mismatch (output: ${output_records}, expected: ${expected_records})")
        result_status[$i]="FAIL"
        result_reason[$i]="bcftools stats — record count mismatch (output: ${output_records}, expected: ${expected_records})"
        result_detail[$i]=""
        result_full_diff[$i]=""
        result_isec_output_only[$i]=""
        result_isec_expected_only[$i]=""
        continue
    fi

    # ── Tier 2: bcftools isec ───────────────────────────────────────────────

    mkdir -p "${isec_dir}"
    if ! bcftools isec -p "${isec_dir}" "${local_output}" "${local_expected}" 2>/dev/null; then
        echo "  [FAIL] ${file} — bcftools isec failed"
        failed=$((failed + 1))
        failures+=("${file}: bcftools isec failed")
        result_status[$i]="FAIL"
        result_reason[$i]="bcftools isec failed"
        result_detail[$i]=""
        result_full_diff[$i]=""
        result_isec_output_only[$i]=""
        result_isec_expected_only[$i]=""
        continue
    fi

    output_only=$(bcftools view -H "${isec_dir}/0000.vcf" 2>/dev/null | wc -l)
    expected_only=$(bcftools view -H "${isec_dir}/0001.vcf" 2>/dev/null | wc -l)
    result_isec_output_only[$i]="${output_only}"
    result_isec_expected_only[$i]="${expected_only}"

    if (( output_only != 0 || expected_only != 0 )); then
        echo "  [FAIL] ${file} — ${output_only} records only in output, ${expected_only} only in expected"
        failed=$((failed + 1))
        failures+=("${file}: ${output_only} output-only, ${expected_only} expected-only")
        result_status[$i]="FAIL"
        result_reason[$i]="bcftools isec — ${output_only} output-only, ${expected_only} expected-only"
        result_detail[$i]=""
        result_full_diff[$i]=""
        continue
    fi

    # ── Tier 3: bcftools query (field-level) ────────────────────────────────

    query_fmt='%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%AD]\n'
    query_output="${TMPDIR}/query_output_${file%.vcf.gz}.txt"
    query_expected="${TMPDIR}/query_expected_${file%.vcf.gz}.txt"
    query_diff="${TMPDIR}/query_diff_${file%.vcf.gz}.txt"

    if ! bcftools query -f "${query_fmt}" "${local_output}" > "${query_output}" 2>/dev/null; then
        echo "  [FAIL] ${file} — bcftools query failed on output"
        failed=$((failed + 1))
        failures+=("${file}: bcftools query failed on output")
        result_status[$i]="FAIL"
        result_reason[$i]="bcftools query failed on output"
        result_detail[$i]=""
        result_full_diff[$i]=""
        continue
    fi

    if ! bcftools query -f "${query_fmt}" "${local_expected}" > "${query_expected}" 2>/dev/null; then
        echo "  [FAIL] ${file} — bcftools query failed on expected"
        failed=$((failed + 1))
        failures+=("${file}: bcftools query failed on expected")
        result_status[$i]="FAIL"
        result_reason[$i]="bcftools query failed on expected"
        result_detail[$i]=""
        result_full_diff[$i]=""
        continue
    fi

    if ! diff "${query_output}" "${query_expected}" > "${query_diff}" 2>/dev/null; then
        diff_lines=$(wc -l < "${query_diff}")
        full_diff=$(cat "${query_diff}")
        truncated_diff=$(head -20 "${query_diff}")

        echo "  [FAIL] ${file} — field-level differences (${diff_lines} diff lines)"
        failed=$((failed + 1))
        failures+=("${file}: field-level diff (${diff_lines} lines)")
        result_status[$i]="FAIL"
        result_reason[$i]="bcftools query — field-level differences (${diff_lines} diff lines)"
        result_detail[$i]="${truncated_diff}"
        result_full_diff[$i]="${full_diff}"
        continue
    fi

    # All three tiers passed
    echo "  [PASS] ${file}"
    passed=$((passed + 1))
    result_status[$i]="PASS"
    result_reason[$i]=""
    result_detail[$i]=""
    result_full_diff[$i]=""
done

# ── Summary ──────────────────────────────────────────────────────────────────

echo ""
echo "Results: ${passed}/${total} passed, ${failed} failed"

if [[ ${#failures[@]} -gt 0 ]]; then
    echo ""
    echo "Failures:"
    for f in "${failures[@]}"; do
        echo "  - ${f}"
    done
fi

# ── Markdown report ─────────────────────────────────────────────────────────

{
    echo "# Integration Test Report"
    echo ""
    echo "**Date:** $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
    echo "**Lambda Function:** ${FUNCTION_NAME}"
    echo "**S3 Bucket:** ${BUCKET}"
    echo "**Input path prefix:** ${INPUT_PREFIX}"
    echo ""
    echo "**Result:** ${passed}/${total} passed, ${failed} failed"
    echo ""
    echo "## Results"
    echo ""
    echo "| File | Input variants | Expected variants | Output variants | Status |"
    echo "|------|---------------:|------------------:|----------------:|--------|"

    for i in "${!FILES[@]}"; do
        file="${FILES[$i]}"
        in_c="${result_input_count[$i]:-—}"
        exp_c="${result_expected_count[$i]:-—}"
        out_c="${result_output_count[$i]:-—}"
        status="${result_status[$i]}"

        if [[ "${status}" == "PASS" ]]; then
            status_cell="PASS"
        else
            status_cell="FAIL — ${result_reason[$i]}"
        fi

        echo "| ${file} | ${in_c} | ${exp_c} | ${out_c} | ${status_cell} |"
    done

    # Failure details section
    has_details=false
    for i in "${!FILES[@]}"; do
        if [[ "${result_status[$i]}" == "FAIL" && -n "${result_reason[$i]}" ]]; then
            has_details=true
            break
        fi
    done

    if [[ "${has_details}" == true ]]; then
        echo ""
        echo "## Failure Details"

        for i in "${!FILES[@]}"; do
            if [[ "${result_status[$i]}" != "FAIL" ]]; then
                continue
            fi

            file="${FILES[$i]}"
            echo ""
            echo "### ${file}"
            echo ""
            echo "**Stage failed:** ${result_reason[$i]}"

            # Show isec result if we got that far
            if [[ -n "${result_isec_output_only[$i]:-}" ]]; then
                echo "**isec result:** PASS (${result_isec_output_only[$i]} output-only, ${result_isec_expected_only[$i]} expected-only)"
            elif [[ "${result_reason[$i]}" == *"isec"* ]]; then
                echo "**isec result:** FAIL"
            fi

            # Show field diff if available
            if [[ -n "${result_detail[$i]}" ]]; then
                echo "**Field differences (first 20 lines):**"
                echo '```'
                echo "${result_detail[$i]}"
                echo '```'

                full="${result_full_diff[$i]}"
                full_lines=$(echo "${full}" | wc -l)
                if (( full_lines > 20 )); then
                    echo "*(truncated — full diff in appendix)*"
                fi
            fi
        done
    fi

    # Appendix for full diffs that exceed 20 lines
    has_appendix=false
    for i in "${!FILES[@]}"; do
        if [[ -n "${result_full_diff[$i]:-}" ]]; then
            full_lines=$(echo "${result_full_diff[$i]}" | wc -l)
            if (( full_lines > 20 )); then
                has_appendix=true
                break
            fi
        fi
    done

    if [[ "${has_appendix}" == true ]]; then
        echo ""
        echo "## Appendix — Full Diffs"

        for i in "${!FILES[@]}"; do
            if [[ -z "${result_full_diff[$i]:-}" ]]; then
                continue
            fi

            full_lines=$(echo "${result_full_diff[$i]}" | wc -l)
            if (( full_lines <= 20 )); then
                continue
            fi

            echo ""
            echo "### ${FILES[$i]}"
            echo ""
            echo '```'
            echo "${result_full_diff[$i]}"
            echo '```'
        done
    fi
} > "${REPORT_FILE}"

echo ""
echo "Report written to ${REPORT_FILE}"

if [[ ${failed} -gt 0 ]]; then
    exit 1
fi
