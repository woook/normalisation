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
