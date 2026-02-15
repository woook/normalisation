terraform {
  required_version = ">= 1.5"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
}

data "aws_caller_identity" "current" {}
data "aws_region" "current" {}

# ---------------------------------------------------------------------------
# ECR Repository
# ---------------------------------------------------------------------------

resource "aws_ecr_repository" "this" {
  name                 = var.project_name
  image_tag_mutability = "MUTABLE"
  force_delete         = true

  image_scanning_configuration {
    scan_on_push = true
  }

  tags = var.tags
}

resource "aws_ecr_lifecycle_policy" "this" {
  repository = aws_ecr_repository.this.name

  policy = jsonencode({
    rules = [{
      rulePriority = 1
      description  = "Keep only the last 5 images"
      selection = {
        tagStatus   = "any"
        countType   = "imageCountMoreThan"
        countNumber = 5
      }
      action = {
        type = "expire"
      }
    }]
  })
}

# ---------------------------------------------------------------------------
# IAM Role for Lambda
# ---------------------------------------------------------------------------

resource "aws_iam_role" "lambda" {
  name = "${var.project_name}-lambda"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Action = "sts:AssumeRole"
      Effect = "Allow"
      Principal = {
        Service = "lambda.amazonaws.com"
      }
    }]
  })

  tags = var.tags
}

resource "aws_iam_role_policy_attachment" "lambda_basic" {
  role       = aws_iam_role.lambda.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole"
}

resource "aws_iam_role_policy" "lambda_s3" {
  name = "${var.project_name}-s3-access"
  role = aws_iam_role.lambda.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid    = "ReadInputBucket"
        Effect = "Allow"
        Action = [
          "s3:GetObject",
        ]
        Resource = "${data.aws_s3_bucket.input.arn}/${var.input_prefix}*"
      },
      {
        Sid    = "WriteOutputBucket"
        Effect = "Allow"
        Action = [
          "s3:PutObject",
        ]
        Resource = "${data.aws_s3_bucket.input.arn}/${var.output_prefix}*"
      },
      {
        Sid    = "ReadGenomeRef"
        Effect = "Allow"
        Action = [
          "s3:GetObject",
        ]
        Resource = [
          "arn:aws:s3:::${var.genome_ref_bucket}/${var.genome_ref_key}",
          "arn:aws:s3:::${var.genome_ref_bucket}/${var.genome_ref_key}.fai",
          "arn:aws:s3:::${var.genome_ref_bucket}/${var.genome_ref_key}.gzi",
        ]
      },
    ]
  })
}

# ---------------------------------------------------------------------------
# S3 Bucket (data source — bucket must already exist)
# ---------------------------------------------------------------------------

data "aws_s3_bucket" "input" {
  bucket = var.input_bucket_name
}

# ---------------------------------------------------------------------------
# Lambda Function
# ---------------------------------------------------------------------------

resource "aws_lambda_function" "normalise" {
  function_name = var.project_name
  role          = aws_iam_role.lambda.arn
  package_type  = "Image"
  image_uri     = "${aws_ecr_repository.this.repository_url}:${var.ecr_image_tag}"
  timeout       = var.lambda_timeout
  memory_size   = var.lambda_memory_mb

  ephemeral_storage {
    size = var.lambda_ephemeral_storage_mb
  }

  environment {
    variables = {
      GENOME_REF_BUCKET = var.genome_ref_bucket
      GENOME_REF_KEY    = var.genome_ref_key
      OUTPUT_PREFIX     = var.output_prefix
    }
  }

  tags = var.tags
}

# ---------------------------------------------------------------------------
# S3 Event Notification → Lambda
# ---------------------------------------------------------------------------

resource "aws_lambda_permission" "s3" {
  statement_id  = "AllowS3Invoke"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.normalise.function_name
  principal     = "s3.amazonaws.com"
  source_arn    = data.aws_s3_bucket.input.arn
  source_account = data.aws_caller_identity.current.account_id
}

resource "aws_s3_bucket_notification" "input" {
  bucket = data.aws_s3_bucket.input.id

  lambda_function {
    lambda_function_arn = aws_lambda_function.normalise.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = var.input_prefix
    filter_suffix       = ".vcf.gz"
  }

  depends_on = [aws_lambda_permission.s3]
}
