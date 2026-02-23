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
input_bucket_name = "my-vcf-data"
genome_ref_bucket = "my-reference-genomes"
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

Upload a VCF to the `input/` prefix in your bucket. Both gzipped (`.vcf.gz`) and uncompressed (`.vcf`) inputs are accepted. The Lambda triggers automatically:

```bash
aws s3 cp sample.vcf.gz s3://my-vcf-data/input/sample.vcf.gz
# or
aws s3 cp sample.vcf s3://my-vcf-data/input/sample.vcf
```

The normalised file appears at the `output/` prefix. Output is always bgzipped (`.vcf.gz`), regardless of whether the input was compressed:

```bash
# Check it arrived
aws s3 ls s3://my-vcf-data/output/sample.vcf.gz

# Download it
aws s3 cp s3://my-vcf-data/output/sample.vcf.gz normalised_sample.vcf.gz
```

### Manual — re-process a file

Use the helper script to re-invoke the Lambda for a specific file:

```bash
./scripts/invoke.sh my-vcf-data input/sample.vcf.gz
```

Or invoke directly with the AWS CLI:

```bash
aws lambda invoke \
  --function-name "$FUNCTION_NAME" \
  --payload '{"bucket": "my-vcf-data", "key": "input/sample.vcf.gz"}' \
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
| `extra_s3_prefixes` | Extra read/write S3 prefix pairs for Lambda (e.g. testing) | `[]` |

## Testing

### Unit tests

```bash
pytest tests/
```

### Integration tests

The integration test script invokes the Lambda on all test VCFs stored in S3 and compares outputs against expected files using three levels of comparison:

1. **`bcftools stats`** — record count sanity check (catches obvious mismatches early)
2. **`bcftools isec`** — site-level comparison (CHROM/POS/REF/ALT identity)
3. **`bcftools query` + `diff`** — field-level comparison (GT, DP, AD values — catches e.g. incorrect AD splits after multiallelic decomposition)

A file passes only if all three tiers pass. After the run, a markdown report is written to `integration_report.md` (configurable via `REPORT_FILE`) with a results table, failure details, and full diffs in an appendix.

Test data lives under separate S3 prefixes (`test/input/`, `test/expected/`), completely separate from the production `input/` → `output/` flow.

#### Local dependencies

The integration test script runs `bcftools` locally to compare output and expected files. You need:

- **bcftools** (>= 1.13) — for `bcftools stats`, `bcftools isec`, `bcftools query`, `bcftools index`, and `bcftools view`
- **bgzip** (from htslib) — bundled with most bcftools installations
- **jq** — used to construct the Lambda invocation payload
- **AWS CLI v2** — for S3 downloads and Lambda invocation
- **diff** — standard Unix diff (coreutils)

On Ubuntu/Debian:

```bash
sudo apt install bcftools jq
```

On macOS (Homebrew):

```bash
brew install bcftools jq
```

AWS CLI v2 must be installed separately — see the [AWS CLI installation guide](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html). `diff` is preinstalled on all standard Linux and macOS systems.

#### Setup

1. Grant the Lambda access to the test prefixes by adding the following to `terraform.tfvars` and running `terraform apply`:

```hcl
extra_s3_prefixes = [
  {
    read_prefix  = "test/input/"
    write_prefix = "test/output/"
  }
]
```

2. Upload test inputs and expected outputs to S3:

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
| `REPORT_FILE` | Env var | `integration_report.md` |

Example:

```bash
MAX_PARALLEL=20 ./scripts/integration_test.sh my-vcf-data
```

## AWS permissions

### Deployment

The user or CI role running `terraform apply` and pushing the container image needs:

| Service | Permissions | Reason |
|---------|-------------|--------|
| ECR | `ecr:CreateRepository`, `ecr:DeleteRepository`, `ecr:PutLifecyclePolicy`, `ecr:DescribeRepositories`, `ecr:ListTagsForResource`, `ecr:TagResource` | Create and manage the container repository |
| ECR (image push) | `ecr:GetAuthorizationToken`, `ecr:BatchCheckLayerAvailability`, `ecr:PutImage`, `ecr:InitiateLayerUpload`, `ecr:UploadLayerPart`, `ecr:CompleteLayerUpload` | Authenticate Docker and push images |
| Lambda | `lambda:CreateFunction`, `lambda:UpdateFunctionCode`, `lambda:UpdateFunctionConfiguration`, `lambda:DeleteFunction`, `lambda:GetFunction`, `lambda:AddPermission`, `lambda:RemovePermission`, `lambda:TagResource`, `lambda:ListTags` | Create and update the Lambda function |
| IAM | `iam:CreateRole`, `iam:DeleteRole`, `iam:AttachRolePolicy`, `iam:DetachRolePolicy`, `iam:PutRolePolicy`, `iam:DeleteRolePolicy`, `iam:GetRole`, `iam:GetRolePolicy`, `iam:PassRole`, `iam:ListRolePolicies`, `iam:ListAttachedRolePolicies`, `iam:ListInstanceProfilesForRole`, `iam:TagRole` | Manage the Lambda execution role |
| S3 | `s3:GetBucketNotification`, `s3:PutBucketNotification` | Configure the S3 event trigger |
| S3 (data source) | `s3:ListBucket`, `s3:GetBucketLocation` | Terraform data source to reference the existing bucket |
| STS | `sts:GetCallerIdentity` | Terraform uses this to determine the account ID |

### Runtime (day-to-day use)

Users who upload VCFs or manually invoke the Lambda need:

| Service | Permissions | Reason |
|---------|-------------|--------|
| S3 | `s3:PutObject` on `input/*` | Upload input VCFs |
| S3 | `s3:GetObject` on `output/*` | Download normalised results |
| S3 | `s3:ListBucket` | List objects in input/output prefixes |
| Lambda | `lambda:InvokeFunction` | Manual invocation via `invoke.sh` or AWS CLI |
| CloudWatch Logs | `logs:FilterLogEvents`, `logs:GetLogEvents` | View Lambda logs for monitoring |

### Integration testing

In addition to the runtime permissions above, the integration test script needs:

| Service | Permissions | Reason |
|---------|-------------|--------|
| S3 | `s3:GetObject` on `test/input/*`, `test/expected/*`, `test/output/*` | Download test inputs, expected files, and Lambda outputs |
| S3 | `s3:DeleteObject` on `test/output/*` | Clean previous test outputs before each run |
| S3 | `s3:ListBucket` (with prefix `test/input/`) | Discover test VCF files |
| Lambda | `lambda:InvokeFunction` (async) | Invoke the Lambda for each test file |

## Normalisation command

The pipeline runs:

```bash
bcftools norm -Oz -f genome.fa -m -any --keep-sum AD input.vcf.gz -o output.vcf.gz
```

- `-m -any` — split multiallelic sites into biallelic records
- `--keep-sum AD` — maintain the allelic depth sum when splitting
- `-f genome.fa` — left-align and normalise indels against the reference

## Further reading

See [technical_walkthrough.md](technical_walkthrough.md) for a detailed walkthrough of the codebase including the handler source, Dockerfile, Terraform resources, and test output.
