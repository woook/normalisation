"""Unit tests for the VCF normalisation Lambda handler."""

import json
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Set required env vars before importing the handler
os.environ.setdefault("GENOME_REF_BUCKET", "test-genome-bucket")
os.environ.setdefault("GENOME_REF_KEY", "genomes/hg38/genome.fa")
os.environ.setdefault("OUTPUT_PREFIX", "output/")

from handler import (
    _download_genome,
    _parse_event,
    _run_bcftools_norm,
    _upload_output,
    lambda_handler,
)


# ---------------------------------------------------------------------------
# Event parsing
# ---------------------------------------------------------------------------


class TestParseEvent:
    """Tests for _parse_event handling of S3 and manual payloads."""

    def test_s3_event(self):
        event = {
            "Records": [
                {
                    "s3": {
                        "bucket": {"name": "my-bucket"},
                        "object": {"key": "input/sample.vcf.gz"},
                    }
                }
            ]
        }
        bucket, key = _parse_event(event)
        assert bucket == "my-bucket"
        assert key == "input/sample.vcf.gz"

    def test_s3_event_url_encoded_key(self):
        event = {
            "Records": [
                {
                    "s3": {
                        "bucket": {"name": "my-bucket"},
                        "object": {"key": "input/my+sample.vcf.gz"},
                    }
                }
            ]
        }
        bucket, key = _parse_event(event)
        assert key == "input/my sample.vcf.gz"

    def test_manual_event(self):
        event = {"bucket": "my-bucket", "key": "input/sample.vcf.gz"}
        bucket, key = _parse_event(event)
        assert bucket == "my-bucket"
        assert key == "input/sample.vcf.gz"

    def test_missing_fields_raises(self):
        with pytest.raises((KeyError, IndexError)):
            _parse_event({})


# ---------------------------------------------------------------------------
# Genome download
# ---------------------------------------------------------------------------


class TestDownloadGenome:
    """Tests for _download_genome with uncompressed and bgzipped genomes."""

    @patch("handler.s3")
    def test_uncompressed_genome(self, mock_s3, tmp_path):
        with patch("handler.WORK_DIR", tmp_path), \
             patch("handler.GENOME_REF_BUCKET", "ref-bucket"), \
             patch("handler.GENOME_REF_KEY", "genomes/genome.fa"):
            result = _download_genome()

        assert result == tmp_path / "genome.fa"
        assert mock_s3.download_file.call_count == 2
        mock_s3.download_file.assert_any_call("ref-bucket", "genomes/genome.fa", str(tmp_path / "genome.fa"))
        mock_s3.download_file.assert_any_call("ref-bucket", "genomes/genome.fa.fai", str(tmp_path / "genome.fa.fai"))

    @patch("handler.s3")
    def test_bgzipped_genome_downloads_gzi(self, mock_s3, tmp_path):
        with patch("handler.WORK_DIR", tmp_path), \
             patch("handler.GENOME_REF_BUCKET", "ref-bucket"), \
             patch("handler.GENOME_REF_KEY", "genomes/genome.fa.gz"):
            result = _download_genome()

        assert result == tmp_path / "genome.fa.gz"
        assert mock_s3.download_file.call_count == 3
        mock_s3.download_file.assert_any_call("ref-bucket", "genomes/genome.fa.gz", str(tmp_path / "genome.fa.gz"))
        mock_s3.download_file.assert_any_call("ref-bucket", "genomes/genome.fa.gz.fai", str(tmp_path / "genome.fa.gz.fai"))
        mock_s3.download_file.assert_any_call("ref-bucket", "genomes/genome.fa.gz.gzi", str(tmp_path / "genome.fa.gz.gzi"))


# ---------------------------------------------------------------------------
# bcftools norm subprocess
# ---------------------------------------------------------------------------


class TestBcftoolsNorm:
    """Tests for _run_bcftools_norm subprocess execution."""

    @patch("handler.subprocess.run")
    def test_success(self, mock_run, tmp_path):
        input_path = tmp_path / "sample.vcf.gz"
        input_path.touch()
        genome_path = tmp_path / "genome.fa"
        genome_path.touch()

        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="Lines total/split/joined: 100/10/5\n")

        # Patch WORK_DIR so output lands in tmp_path
        with patch("handler.WORK_DIR", tmp_path):
            output = _run_bcftools_norm(input_path, genome_path)

        assert output == tmp_path / f"normalised_{input_path.name}"
        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "bcftools"
        assert "-f" in cmd
        assert "--keep-sum" in cmd

    @patch("handler.subprocess.run")
    def test_failure_raises(self, mock_run, tmp_path):
        input_path = tmp_path / "sample.vcf.gz"
        input_path.touch()
        genome_path = tmp_path / "genome.fa"
        genome_path.touch()

        mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="error: genome mismatch")

        with patch("handler.WORK_DIR", tmp_path):
            with pytest.raises(RuntimeError, match="bcftools norm failed"):
                _run_bcftools_norm(input_path, genome_path)


# ---------------------------------------------------------------------------
# Upload output
# ---------------------------------------------------------------------------


class TestUploadOutput:
    """Tests for _upload_output S3 key derivation logic."""

    @patch("handler.s3")
    def test_standard_input_prefix(self, mock_s3, tmp_path):
        """input/sample.vcf.gz should produce output/sample.vcf.gz."""
        output_path = tmp_path / "normalised_sample.vcf.gz"
        output_path.touch()

        key = _upload_output("my-bucket", "input/sample.vcf.gz", output_path)
        assert key == "output/sample.vcf.gz"
        mock_s3.upload_file.assert_called_once_with(
            str(output_path), "my-bucket", "output/sample.vcf.gz"
        )

    @patch("handler.s3")
    def test_nested_input_prefix(self, mock_s3, tmp_path):
        """test/input/sample.vcf.gz should produce test/output/sample.vcf.gz."""
        output_path = tmp_path / "normalised_sample.vcf.gz"
        output_path.touch()

        key = _upload_output("my-bucket", "test/input/sample.vcf.gz", output_path)
        assert key == "test/output/sample.vcf.gz"
        mock_s3.upload_file.assert_called_once_with(
            str(output_path), "my-bucket", "test/output/sample.vcf.gz"
        )

    @patch("handler.s3")
    def test_multiple_input_segments_replaces_last(self, mock_s3, tmp_path):
        """When multiple /input/ segments exist, only the last is replaced."""
        output_path = tmp_path / "normalised_sample.vcf.gz"
        output_path.touch()

        key = _upload_output("my-bucket", "input/subdir/input/sample.vcf.gz", output_path)
        assert key == "input/subdir/output/sample.vcf.gz"
        mock_s3.upload_file.assert_called_once_with(
            str(output_path), "my-bucket", "input/subdir/output/sample.vcf.gz"
        )

    @patch("handler.s3")
    def test_no_input_segment_uses_output_prefix(self, mock_s3, tmp_path):
        """Keys without /input/ fall back to OUTPUT_PREFIX."""
        output_path = tmp_path / "normalised_sample.vcf.gz"
        output_path.touch()

        key = _upload_output("my-bucket", "uploads/sample.vcf.gz", output_path)
        assert key == "output/sample.vcf.gz"
        mock_s3.upload_file.assert_called_once_with(
            str(output_path), "my-bucket", "output/sample.vcf.gz"
        )


# ---------------------------------------------------------------------------
# Full handler (integration-style with mocks)
# ---------------------------------------------------------------------------


class TestLambdaHandler:
    """End-to-end handler tests with mocked dependencies."""

    @patch("handler._cleanup")
    @patch("handler._upload_output", return_value="output/sample.vcf.gz")
    @patch("handler._run_bcftools_norm")
    @patch("handler._download_genome")
    @patch("handler._download_input")
    @patch("handler._setup_work_dir")
    def test_full_flow(
        self,
        mock_setup,
        mock_dl_input,
        mock_dl_genome,
        mock_norm,
        mock_upload,
        mock_cleanup,
    ):
        event = {"bucket": "my-bucket", "key": "input/sample.vcf.gz"}
        result = lambda_handler(event, None)

        assert result["statusCode"] == 200
        body = json.loads(result["body"])
        assert body["input"] == "s3://my-bucket/input/sample.vcf.gz"
        assert body["output"] == "s3://my-bucket/output/sample.vcf.gz"

        mock_setup.assert_called_once()
        mock_dl_input.assert_called_once_with("my-bucket", "input/sample.vcf.gz")
        mock_dl_genome.assert_called_once()
        mock_norm.assert_called_once()
        mock_upload.assert_called_once()
        mock_cleanup.assert_called_once()

    @patch("handler._cleanup")
    @patch("handler._download_input", side_effect=Exception("S3 error"))
    @patch("handler._setup_work_dir")
    def test_cleanup_on_error(self, mock_setup, mock_dl, mock_cleanup):
        event = {"bucket": "my-bucket", "key": "input/sample.vcf.gz"}

        with pytest.raises(Exception, match="S3 error"):
            lambda_handler(event, None)

        # Cleanup should still be called via finally
        mock_cleanup.assert_called_once()
