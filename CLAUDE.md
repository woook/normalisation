# VCF Normalisation Pipeline

## Quick reference

- **Language**: Python 3.12 (Lambda handler), HCL (Terraform)
- **Test command**: `pytest tests/`
- **Terraform validate**: `cd terraform && terraform init && terraform validate`
- **Docker build**: `docker build -t vcf-normalisation .`

## Architecture

S3 event (input/*.vcf.gz) → Lambda (container image with bcftools) → S3 (output/)

The Lambda downloads the input VCF and a per-account reference genome from S3,
runs `bcftools norm -Oz -f genome.fa -m -any --keep-sum AD`, and uploads the result.

## Key files

- `lambda/handler.py` — Lambda entry point
- `Dockerfile` — Amazon Linux + bcftools/htslib container image
- `terraform/main.tf` — Infrastructure (ECR, Lambda, S3 notification, IAM)
- `scripts/invoke.sh` — Manual invocation helper

## Conventions

- The same Terraform module is deployed to 7 AWS accounts with different tfvars
- Input bucket must already exist; Terraform uses a data source to reference it
- Output goes to `output/` prefix in the same bucket by default
