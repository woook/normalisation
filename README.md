# VCF Normalisation Pipeline

Serverless pipeline that normalises VCF files using `bcftools norm`. Deployed as a Lambda container image triggered by S3 uploads.

## Architecture

```text
S3 (input/)  →  S3 Event Notification  →  Lambda (bcftools container)  →  S3 (output/)
```

Each of the 7 groups deploys the same Terraform module into their own AWS account.

## Prerequisites

- AWS CLI v2 configured with credentials for the target account
- Terraform >= 1.5
- Docker
- An existing S3 bucket for VCF uploads (the "input bucket")
- A reference genome uploaded to an S3 bucket — either uncompressed (`.fa` + `.fa.fai`) or bgzipped (`.fa.gz` + `.fa.gz.fai` + `.fa.gz.gzi`)

## Deployment

Deployment is a two-pass process: Terraform creates the ECR repository first, then you push the container image and apply again to update the Lambda.

Before starting, ensure your AWS credentials are valid and not expired:

```bash
aws sts get-caller-identity
```

If using SSO, run `aws sso login` first. Terraform and the AWS CLI both require active credentials for every step below.

### Step 1 — Configure Terraform variables

```bash
cd terraform
cp terraform.tfvars.example terraform.tfvars
```

Edit `terraform.tfvars` with your account-specific values:

```hcl
input_bucket_name = "my-group-vcf-data"
genome_ref_bucket = "my-group-reference-genomes"
genome_ref_key    = "genomes/hg38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
```

The following index files must exist alongside the genome in the same bucket:

- **Uncompressed genome** (`.fa`): requires `.fa.fai`
- **Bgzipped genome** (`.fa.gz`): requires `.fa.gz.fai` and `.fa.gz.gzi`

The pipeline detects the format from the file extension and downloads the appropriate indices automatically.

### Step 2 — Create the ECR repository

```bash
terraform init
terraform apply -target=aws_ecr_repository.this -target=aws_ecr_lifecycle_policy.this
```

This creates the ECR repository without attempting to create the Lambda (which needs the image to exist first).

Grab the ECR repository URL from the output:

```bash
ECR_REPO=$(terraform output -raw ecr_repository_url)
```

### Step 3 — Build and push the container image

> **Note:** Docker commands require root privileges. Prefix with `sudo` or run as root.

```bash
# From the project root
cd ..

# Build the image
sudo docker build -t vcf-normalisation .

# Authenticate Docker to ECR
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
REGION=$(aws configure get region)
aws ecr get-login-password --region "$REGION" \
  | sudo docker login --username AWS --password-stdin "$ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com"

# Tag and push
sudo docker tag vcf-normalisation:latest "$ECR_REPO:latest"
sudo docker push "$ECR_REPO:latest"
```

### Step 4 — Full Terraform apply (deploys Lambda)

```bash
cd terraform
terraform apply
```

This creates the Lambda function, IAM roles, and S3 event notification using the image you just pushed.

### Updating the Lambda

When you update the handler code or bcftools version:

```bash
sudo docker build -t vcf-normalisation .
sudo docker tag vcf-normalisation:latest "$ECR_REPO:latest"
sudo docker push "$ECR_REPO:latest"

# Force Lambda to pick up the new image
FUNCTION_NAME=$(cd terraform && terraform output -raw lambda_function_name)
aws lambda update-function-code \
  --function-name "$FUNCTION_NAME" \
  --image-uri "$ECR_REPO:latest"
```

### Tearing down

```bash
cd terraform
terraform destroy
```

## Running the normalisation

### Automatic — upload a file

Upload a compressed VCF to the `input/` prefix in your bucket. The Lambda triggers automatically:

```bash
aws s3 cp sample.vcf.gz s3://my-group-vcf-data/input/sample.vcf.gz
```

The normalised file appears at the `output/` prefix:

```bash
# Check it arrived
aws s3 ls s3://my-group-vcf-data/output/sample.vcf.gz

# Download it
aws s3 cp s3://my-group-vcf-data/output/sample.vcf.gz normalised_sample.vcf.gz
```

### Manual — re-process a file

Use the helper script to re-invoke the Lambda for a specific file:

```bash
./scripts/invoke.sh my-group-vcf-data input/sample.vcf.gz
```

Or invoke directly with the AWS CLI:

```bash
aws lambda invoke \
  --function-name "$FUNCTION_NAME" \
  --payload '{"bucket": "my-group-vcf-data", "key": "input/sample.vcf.gz"}' \
  --cli-binary-format raw-in-base64-out \
  /dev/stdout
```

### Monitoring

Set the function name from Terraform output (or use the `FUNCTION_NAME` env var from earlier steps):

```bash
FUNCTION_NAME=$(cd terraform && terraform output -raw lambda_function_name)
```

View Lambda logs in CloudWatch:

```bash
aws logs tail "/aws/lambda/$FUNCTION_NAME" --follow
```

Check for invocation errors:

```bash
aws lambda get-function --function-name "$FUNCTION_NAME" \
  --query 'Configuration.{State:State,LastModified:LastModified}'
```

## Configuration

See `terraform/variables.tf` for all options. Key variables:

| Variable | Description | Default |
|---|---|---|
| `input_bucket_name` | S3 bucket for VCF uploads | (required) |
| `genome_ref_bucket` | S3 bucket with reference genome | (required) |
| `genome_ref_key` | S3 key for the genome (`.fa` or `.fa.gz`) | (required) |
| `input_prefix` | S3 prefix that triggers Lambda | `input/` |
| `output_prefix` | S3 prefix for normalised output | `output/` |
| `lambda_memory_mb` | Lambda memory (MB) | 2048 |
| `lambda_timeout` | Lambda timeout (seconds) | 600 |
| `lambda_ephemeral_storage_mb` | Ephemeral `/tmp` storage (MB) | 4096 |
| `ecr_image_tag` | Container image tag to deploy | `latest` |

## Testing

### Unit tests

```bash
pytest tests/
```

### Integration tests

The integration test script invokes the Lambda on all test VCFs stored in S3 and compares outputs against pre-normalised expected files using `bcftools isec`.

Test data lives under separate S3 prefixes (`test/input/`, `test/expected/`), completely separate from the production `input/` → `output/` flow.

#### Setup

Upload test inputs and expected outputs to S3:

```bash
aws s3 sync ./test_vcfs/ s3://my-vcf-data/test/input/
aws s3 sync ./expected_vcfs/ s3://my-vcf-data/test/expected/
```

#### Running

```bash
./scripts/integration_test.sh <bucket> [input_prefix] [expected_prefix]
```

| Parameter | Source | Default |
|---|---|---|
| `BUCKET` | Arg 1 | (required) |
| `INPUT_PREFIX` | Arg 2 | `test/input/` |
| `EXPECTED_PREFIX` | Arg 3 | `test/expected/` |
| `FUNCTION_NAME` | Env var | `vcf-normalisation` |
| `MAX_PARALLEL` | Env var | `10` |
| `POLL_TIMEOUT` | Env var | `120` (seconds per file) |

Example:

```bash
FUNCTION_NAME=my-vcf-lambda MAX_PARALLEL=20 ./scripts/integration_test.sh my-vcf-data
```

## Normalisation command

The pipeline runs:

```bash
bcftools norm -f genome.fa -m -any --keep-sum AD input.vcf.gz -o output.vcf.gz
```

- `-m -any` — split multiallelic sites into biallelic records
- `--keep-sum AD` — maintain the allelic depth sum when splitting
- `-f genome.fa` — left-align and normalise indels against the reference
