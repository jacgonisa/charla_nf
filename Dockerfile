# Stage 1: The Rust Builder Stage
FROM rust:latest AS rust-builder

WORKDIR /usr/src/app

COPY rust/combine_segment_fast_indixes /usr/src/combine_segment_fast_indixes
COPY rust/kmer_processor /usr/src/kmer_processor

RUN cargo install --path /usr/src/combine_segment_fast_indixes --root /usr/local
RUN cargo install --path /usr/src/kmer_processor --root /usr/local

#-----------------------------------------------------------------------------------------------------------------------------------

# Stage 2: The KMC Builder Stage
FROM ubuntu:24.04 AS kmc-builder

# Install necessary tools to add the deadsnakes PPA.
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    && rm -rf /var/lib/apt/lists/*

# Add the deadsnakes PPA for Python 3.10.
RUN add-apt-repository ppa:deadsnakes/ppa

# Install all build dependencies.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    git \
    python3.10 \
    python3.10-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Create a symlink for python3 to point to python3.10.
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 1

# Copy and build the KMC source code.
COPY kmc /usr/src/kmc
WORKDIR /usr/src/kmc
RUN make -j8

#-----------------------------------------------------------------------------------------------------------------------------------

# Stage 3: The Final Pipeline Image
FROM continuumio/miniconda3:latest

WORKDIR /app

# Copy the Conda environment file and create the environment.
COPY containers/environment_charla.yml .
RUN conda env create -f environment_charla.yml

# Add the new Conda environment and executables to the system's PATH.
ENV PATH="/opt/conda/envs/charla_nf/bin:$PATH"

# Copy the compiled Rust executables from the 'rust-builder' stage.
COPY --from=rust-builder /usr/local/bin/hybrid_profile_toolkit /usr/local/bin/
COPY --from=rust-builder /usr/local/bin/kmer_processor /usr/local/bin/

# Copy the compiled KMC executables from the 'kmc-builder' stage.
COPY --from=kmc-builder /usr/src/kmc/bin/kmc /usr/local/bin/
COPY --from=kmc-builder /usr/src/kmc/bin/kmc_tools /usr/local/bin/

WORKDIR /data
