# VCF Normalisation Pipeline

*2026-02-19T22:40:27Z by Showboat 0.6.0*
<!-- showboat-id: ee047990-1b34-4542-b149-b846ac76e849 -->

A serverless pipeline that normalises VCF files using `bcftools norm`. An S3 upload to the `input/` prefix triggers a Lambda function (packaged as a container image with bcftools), which left-aligns indels, splits multiallelic sites, and preserves allelic depth sums. Both gzipped (`.vcf.gz`) and uncompressed (`.vcf`) inputs are accepted; output is always bgzipped. The normalised VCF is written to the `output/` prefix in the same bucket. The same Terraform module is deployed across 7 AWS accounts.

## Project Structure

```bash
find . -type f \( -name '*.py' -o -name '*.tf' -o -name '*.sh' -o -name 'Dockerfile' \) | grep -v __pycache__ | grep -v .terraform | sort
```

```output
./Dockerfile
./lambda/handler.py
./scripts/integration_test.sh
./scripts/invoke.sh
./tests/conftest.py
./tests/test_handler.py
```

```bash
ls terraform/*.tf
```

```output
terraform/main.tf
terraform/outputs.tf
terraform/variables.tf
```

## Lambda Handler

The core normalisation logic lives in `lambda/handler.py`. It downloads the input VCF and reference genome from S3, runs `bcftools norm`, and uploads the result. The output path is derived by replacing the last `/input/` segment with `/output/`.

<details>
<summary><code>cat lambda/handler.py</code></summary>

```bash
cat lambda/handler.py
```

```output
"""Lambda handler for VCF normalisation using bcftools."""

import json
import logging
import os
import shutil
import subprocess
import urllib.parse
from pathlib import Path

import boto3

logger = logging.getLogger()
logger.setLevel(logging.INFO)

s3 = boto3.client("s3")

WORK_DIR = Path("/tmp/vcf_work")
GENOME_REF_BUCKET = os.environ["GENOME_REF_BUCKET"]
GENOME_REF_KEY = os.environ["GENOME_REF_KEY"]
OUTPUT_PREFIX = os.environ.get("OUTPUT_PREFIX", "output/")


def lambda_handler(event, context):
    """Entry point for both S3 event triggers and manual invocations.

    Manual invocation payload:
        {"bucket": "my-bucket", "key": "input/sample.vcf.gz"}
    """
    try:
        bucket, key = _parse_event(event)
    except (KeyError, IndexError) as exc:
        logger.error("Failed to parse event: %s", exc)
        raise

    logger.info("Processing s3://%s/%s", bucket, key)

    _setup_work_dir()

    try:
        input_path = _download_input(bucket, key)
        genome_path = _download_genome()
        output_path = _run_bcftools_norm(input_path, genome_path)
        output_key = _upload_output(bucket, key, output_path)
    finally:
        _cleanup()

    logger.info("Normalised file uploaded to s3://%s/%s", bucket, output_key)
    return {
        "statusCode": 200,
        "body": json.dumps(
            {"input": f"s3://{bucket}/{key}", "output": f"s3://{bucket}/{output_key}"}
        ),
    }


def _parse_event(event):
    """Extract bucket and key from S3 event or manual payload."""
    if "Records" in event:
        record = event["Records"][0]["s3"]
        bucket = record["bucket"]["name"]
        key = urllib.parse.unquote_plus(record["object"]["key"])
    else:
        bucket = event["bucket"]
        key = event["key"]
    return bucket, key


def _setup_work_dir():
    """Create the temporary working directory for downloads and processing."""
    WORK_DIR.mkdir(parents=True, exist_ok=True)


def _download_input(bucket, key):
    """Download the input VCF file from S3."""
    filename = Path(key).name
    local_path = WORK_DIR / filename
    logger.info("Downloading input: s3://%s/%s -> %s", bucket, key, local_path)
    s3.download_file(bucket, key, str(local_path))
    return local_path


def _download_genome():
    """Download the reference genome, .fai index, and .gzi index (if bgzipped) from S3."""
    genome_filename = Path(GENOME_REF_KEY).name
    genome_local = WORK_DIR / genome_filename
    fai_key = GENOME_REF_KEY + ".fai"
    fai_local = WORK_DIR / (genome_filename + ".fai")

    logger.info(
        "Downloading genome: s3://%s/%s -> %s",
        GENOME_REF_BUCKET,
        GENOME_REF_KEY,
        genome_local,
    )
    s3.download_file(GENOME_REF_BUCKET, GENOME_REF_KEY, str(genome_local))

    logger.info(
        "Downloading genome index: s3://%s/%s -> %s",
        GENOME_REF_BUCKET,
        fai_key,
        fai_local,
    )
    s3.download_file(GENOME_REF_BUCKET, fai_key, str(fai_local))

    if genome_filename.endswith(".gz"):
        gzi_key = GENOME_REF_KEY + ".gzi"
        gzi_local = WORK_DIR / (genome_filename + ".gzi")
        logger.info(
            "Downloading bgzip index: s3://%s/%s -> %s",
            GENOME_REF_BUCKET,
            gzi_key,
            gzi_local,
        )
        s3.download_file(GENOME_REF_BUCKET, gzi_key, str(gzi_local))

    return genome_local


def _run_bcftools_norm(input_path, genome_path):
    """Run bcftools norm and return the path to the output file.

    Output is always bgzipped (.vcf.gz) regardless of whether the input was
    compressed or uncompressed.
    """
    output_name = input_path.name
    if output_name.endswith(".vcf"):
        output_name = output_name[:-4] + ".vcf.gz"
    output_path = WORK_DIR / f"normalised_{output_name}"

    cmd = [
        "bcftools",
        "norm",
        "-Oz",
        "-f",
        str(genome_path),
        "-m",
        "-any",
        "--keep-sum",
        "AD",
        str(input_path),
        "-o",
        str(output_path),
    ]

    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=540)

    if result.stdout:
        logger.info("bcftools stdout: %s", result.stdout)
    if result.stderr:
        logger.info("bcftools stderr: %s", result.stderr)

    if result.returncode != 0:
        raise RuntimeError(
            f"bcftools norm failed (exit {result.returncode}): {result.stderr}"
        )

    return output_path


def _upload_output(bucket, input_key, output_path):
    """Upload the normalised VCF to S3, deriving the output path from the input key.

    If the input key contains '/input/', the last occurrence is replaced with
    '/output/' to mirror the directory structure. Otherwise falls back to
    OUTPUT_PREFIX + filename.
    """
    if "/input/" in input_key:
        # Replace the last occurrence of /input/ with /output/
        idx = input_key.rfind("/input/")
        output_key = input_key[:idx] + "/output/" + input_key[idx + len("/input/"):]
    else:
        filename = Path(input_key).name
        output_key = f"{OUTPUT_PREFIX}{filename}"

    # Output is always bgzipped; fix the extension if the input was uncompressed
    if output_key.endswith(".vcf"):
        output_key = output_key[:-4] + ".vcf.gz"

    logger.info("Uploading output: %s -> s3://%s/%s", output_path, bucket, output_key)
    s3.upload_file(str(output_path), bucket, output_key)
    return output_key


def _cleanup():
    """Remove all files from the work directory."""
    if WORK_DIR.exists():
        shutil.rmtree(WORK_DIR, ignore_errors=True)
```

</details>

## Container Image

The Lambda runs as a container image built on Amazon Linux with bcftools and htslib compiled from source.

<details>
<summary><code>cat Dockerfile</code></summary>

```bash
cat Dockerfile
```

```output
FROM public.ecr.aws/lambda/python:3.12 AS builder

# Install bcftools and htslib dependencies
RUN dnf install -y \
        autoconf \
        automake \
        bzip2 \
        bzip2-devel \
        gcc \
        libcurl-devel \
        make \
        tar \
        xz-devel \
        zlib-devel \
    && dnf clean all

# Install htslib (provides bgzip, tabix)
ARG HTSLIB_VERSION=1.21
RUN curl -fsSL https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
        | tar -xjf - \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd / && rm -rf htslib-${HTSLIB_VERSION}

# Install bcftools
ARG BCFTOOLS_VERSION=1.21
RUN curl -fsSL https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
        | tar -xjf - \
    && cd bcftools-${BCFTOOLS_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd / && rm -rf bcftools-${BCFTOOLS_VERSION}

# --- Final stage: slim runtime image ---
FROM public.ecr.aws/lambda/python:3.12

# Copy only compiled binaries and shared libraries from builder
COPY --from=builder /usr/local /usr/local
ENV LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

# Copy handler code
COPY lambda/handler.py ${LAMBDA_TASK_ROOT}/handler.py

CMD ["handler.lambda_handler"]
```

</details>

## Normalisation Command

The pipeline runs the following bcftools command:
- `-Oz` — write bgzipped output (output is always `.vcf.gz` regardless of input format)
- `-m -any` — split multiallelic sites into biallelic records
- `--keep-sum AD` — maintain the allelic depth sum when splitting
- `-f genome.fa` — left-align and normalise indels against the reference

```bash
sed -n '/cmd = \[/,/\]/p' lambda/handler.py
```

```output
    cmd = [
        "bcftools",
        "norm",
        "-Oz",
        "-f",
        str(genome_path),
        "-m",
        "-any",
        "--keep-sum",
        "AD",
        str(input_path),
        "-o",
        str(output_path),
    ]
```

## Unit Tests

The test suite validates the handler logic using mocked S3 and subprocess calls.

<details>
<summary><code>pytest tests/ -v</code></summary>

```bash
pytest tests/ -v
```

```output
============================= test session starts ==============================
platform linux -- Python 3.10.19, pytest-9.0.2, pluggy-1.5.0 -- /usr/bin/python3
cachedir: .pytest_cache
rootdir: /home/wook/Documents/normalisation
plugins: anyio-4.6.2.post1
collecting ... collected 17 items

tests/test_handler.py::TestParseEvent::test_s3_event PASSED              [  5%]
tests/test_handler.py::TestParseEvent::test_s3_event_url_encoded_key PASSED [ 11%]
tests/test_handler.py::TestParseEvent::test_manual_event PASSED          [ 17%]
tests/test_handler.py::TestParseEvent::test_missing_fields_raises PASSED [ 23%]
tests/test_handler.py::TestDownloadGenome::test_uncompressed_genome PASSED [ 29%]
tests/test_handler.py::TestDownloadGenome::test_bgzipped_genome_downloads_gzi PASSED [ 35%]
tests/test_handler.py::TestBcftoolsNorm::test_success_gzipped PASSED     [ 41%]
tests/test_handler.py::TestBcftoolsNorm::test_success_uncompressed PASSED [ 47%]
tests/test_handler.py::TestBcftoolsNorm::test_failure_raises PASSED      [ 52%]
tests/test_handler.py::TestUploadOutput::test_standard_input_prefix PASSED [ 58%]
tests/test_handler.py::TestUploadOutput::test_nested_input_prefix PASSED [ 64%]
tests/test_handler.py::TestUploadOutput::test_multiple_input_segments_replaces_last PASSED [ 70%]
tests/test_handler.py::TestUploadOutput::test_no_input_segment_uses_output_prefix PASSED [ 76%]
tests/test_handler.py::TestUploadOutput::test_uncompressed_input_key_becomes_gz PASSED [ 82%]
tests/test_handler.py::TestUploadOutput::test_uncompressed_no_input_segment_becomes_gz PASSED [ 88%]
tests/test_handler.py::TestLambdaHandler::test_full_flow PASSED          [ 94%]
tests/test_handler.py::TestLambdaHandler::test_cleanup_on_error PASSED   [100%]

============================== 17 passed in 0.22s ==============================
```

</details>

## Integration Tests

The integration test script (`scripts/integration_test.sh`) invokes the Lambda on test VCFs in S3 and compares outputs against expected files using three tiers of comparison:

1. **bcftools stats** — record count sanity check
2. **bcftools isec** — site-level (CHROM/POS/REF/ALT) comparison
3. **bcftools query + diff** — field-level (GT/DP/AD) comparison to catch e.g. incorrect AD splits after multiallelic decomposition

A markdown report is generated with a results table, failure details, and full diffs.

```bash
head -21 scripts/integration_test.sh
```

```output
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

```

## Infrastructure

Terraform manages the ECR repository, Lambda function, IAM roles, and S3 event notification. Key variables allow per-account customisation (bucket names, genome reference path, memory/timeout settings).

```bash
grep '^variable' terraform/variables.tf
```

```output
variable "project_name" {
variable "input_bucket_name" {
variable "input_prefix" {
variable "output_prefix" {
variable "genome_ref_bucket" {
variable "genome_ref_key" {
variable "lambda_memory_mb" {
variable "lambda_timeout" {
variable "lambda_ephemeral_storage_mb" {
variable "ecr_image_tag" {
variable "extra_s3_prefixes" {
variable "tags" {
```

```bash
grep '^resource\|^data\|^module' terraform/main.tf
```

```output
data "aws_caller_identity" "current" {}
data "aws_region" "current" {}
resource "aws_ecr_repository" "this" {
resource "aws_ecr_lifecycle_policy" "this" {
resource "aws_iam_role" "lambda" {
resource "aws_iam_role_policy_attachment" "lambda_basic" {
resource "aws_iam_role_policy" "lambda_s3" {
data "aws_s3_bucket" "input" {
resource "aws_lambda_function" "normalise" {
resource "aws_lambda_permission" "s3" {
resource "aws_s3_bucket_notification" "input" {
```

## Input Format Support

The pipeline accepts both gzipped (`.vcf.gz`) and uncompressed (`.vcf`) VCF files. Two changes enable this.

**S3 trigger** — the Terraform notification resource registers two filter blocks, one per suffix, so uploads of either format fire the Lambda:

```bash
grep -A 18 "aws_s3_bucket_notification" terraform/main.tf
```

```output
resource "aws_s3_bucket_notification" "input" {
  bucket = data.aws_s3_bucket.input.id

  lambda_function {
    lambda_function_arn = aws_lambda_function.normalise.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = var.input_prefix
    filter_suffix       = ".vcf.gz"
  }

  lambda_function {
    lambda_function_arn = aws_lambda_function.normalise.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = var.input_prefix
    filter_suffix       = ".vcf"
  }

  depends_on = [aws_lambda_permission.s3]
}
```

**Handler — normalisation and output key** — `_run_bcftools_norm` renames the local output file to `.vcf.gz` when the input was `.vcf`, and passes `-Oz` so bcftools writes bgzipped content. `_upload_output` then ensures the S3 key always ends in `.vcf.gz`:

```bash
grep -A 5 "output_name.endswith" lambda/handler.py
```

```output
    if output_name.endswith(".vcf"):
        output_name = output_name[:-4] + ".vcf.gz"
    output_path = WORK_DIR / f"normalised_{output_name}"

    cmd = [
        "bcftools",
```

```bash
grep -A 2 "Output is always bgzipped" lambda/handler.py
```

```output
    Output is always bgzipped (.vcf.gz) regardless of whether the input was
    compressed or uncompressed.
    """
--
    # Output is always bgzipped; fix the extension if the input was uncompressed
    if output_key.endswith(".vcf"):
        output_key = output_key[:-4] + ".vcf.gz"
```

The test suite covers both input formats end-to-end, verifying that an uncompressed input produces a `.vcf.gz` output path and that the `-Oz` flag is always present in the bcftools command:

```bash
pytest tests/ -v -k "uncompressed"
```

```output
============================= test session starts ==============================
platform linux -- Python 3.10.19, pytest-9.0.2, pluggy-1.5.0 -- /usr/bin/python3
cachedir: .pytest_cache
rootdir: /home/wook/Documents/normalisation
plugins: anyio-4.6.2.post1
collecting ... collected 17 items / 13 deselected / 4 selected

tests/test_handler.py::TestDownloadGenome::test_uncompressed_genome PASSED [ 25%]
tests/test_handler.py::TestBcftoolsNorm::test_success_uncompressed PASSED [ 50%]
tests/test_handler.py::TestUploadOutput::test_uncompressed_input_key_becomes_gz PASSED [ 75%]
tests/test_handler.py::TestUploadOutput::test_uncompressed_no_input_segment_becomes_gz PASSED [100%]

======================= 4 passed, 13 deselected in 0.22s =======================
```
