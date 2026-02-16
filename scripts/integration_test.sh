#!/usr/bin/env bash
#
# Integration test: invoke the VCF normalisation Lambda on test VCFs and
# compare outputs against expected files using bcftools isec.
#
# Usage:
#   ./scripts/integration_test.sh <bucket> [input_prefix] [expected_prefix]
#
# Environment variable overrides:
#   FUNCTION_NAME  Lambda function name        (default: vcf-normalisation)
#   MAX_PARALLEL   Concurrent Lambda invokes   (default: 10)
#   POLL_TIMEOUT   Per-file poll timeout in s   (default: 120)

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

for file in "${FILES[@]}"; do
    # Skip files that never appeared
    if printf '%s\n' "${missing[@]}" 2>/dev/null | grep -qxF "${file}"; then
        echo "  [SKIP] ${file} — output timed out"
        failed=$((failed + 1))
        failures+=("${file}: output timed out")
        continue
    fi

    output_key="${OUTPUT_PREFIX}${file}"
    expected_key="${EXPECTED_PREFIX}${file}"

    local_output="${TMPDIR}/output_${file}"
    local_expected="${TMPDIR}/expected_${file}"
    isec_dir="${TMPDIR}/isec_${file%.vcf.gz}"

    # Download output and expected
    if ! aws s3 cp "s3://${BUCKET}/${output_key}" "${local_output}" > /dev/null 2>&1; then
        echo "  [FAIL] ${file} — could not download output"
        failed=$((failed + 1))
        failures+=("${file}: could not download output")
        continue
    fi

    if ! aws s3 cp "s3://${BUCKET}/${expected_key}" "${local_expected}" > /dev/null 2>&1; then
        echo "  [FAIL] ${file} — could not download expected"
        failed=$((failed + 1))
        failures+=("${file}: could not download expected")
        continue
    fi

    # Index both files for bcftools isec
    bcftools index "${local_output}" 2>/dev/null || true
    bcftools index "${local_expected}" 2>/dev/null || true

    # Run bcftools isec
    mkdir -p "${isec_dir}"
    if ! bcftools isec -p "${isec_dir}" "${local_output}" "${local_expected}" 2>/dev/null; then
        echo "  [FAIL] ${file} — bcftools isec failed"
        failed=$((failed + 1))
        failures+=("${file}: bcftools isec failed")
        continue
    fi

    # Count records unique to each file
    output_only=$(bcftools view -H "${isec_dir}/0000.vcf" 2>/dev/null | wc -l)
    expected_only=$(bcftools view -H "${isec_dir}/0001.vcf" 2>/dev/null | wc -l)

    if (( output_only == 0 && expected_only == 0 )); then
        echo "  [PASS] ${file}"
        passed=$((passed + 1))
    else
        echo "  [FAIL] ${file} — ${output_only} records only in output, ${expected_only} only in expected"
        failed=$((failed + 1))
        failures+=("${file}: ${output_only} output-only, ${expected_only} expected-only")
    fi
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
    exit 1
fi
