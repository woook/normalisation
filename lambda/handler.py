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
