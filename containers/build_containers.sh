#!/bin/bash
# build_containers.sh - Build Docker and convert to Singularity

set -e

# Configuration
IMAGE_NAME="charla-nf"
VERSION="1.0.0"
DOCKER_TAG="${IMAGE_NAME}:${VERSION}"
DOCKER_LATEST="${IMAGE_NAME}:latest"
SINGULARITY_IMAGE="${IMAGE_NAME}_${VERSION}.sif"

echo "========================================"
echo "Building CHARLA-NF Container Images"
echo "========================================"

# Check if Docker is available
if ! command -v docker &> /dev/null; then
    echo "ERROR: Docker is not installed or not in PATH"
    exit 1
fi

# Build Docker image
echo ""
echo "Step 1: Building Docker image..."
echo "----------------------------------------"
docker build \
    --tag ${DOCKER_TAG} \
    --tag ${DOCKER_LATEST} \
    --file charla_nf/Dockerfile \
    .

echo ""
echo "✅ Docker image built successfully: ${DOCKER_TAG}"

# Test Docker image
echo ""
echo "Step 2: Testing Docker image..."
echo "----------------------------------------"
sudo docker run --rm ${DOCKER_LATEST} bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate environment_charla && \                                      
    echo '=== Testing installations ===' && \
    hybrid_profile_toolkit --version || echo 'hybrid_profile_toolkit: OK' && \
    kmer_processor --version || echo 'kmer_processor: OK' && \
    kmc || echo 'KMC: OK' && \
    seqkit version && \
    minimap2 --version && \
    samtools --version && \
    R --version && \
    python --version && \
    echo '=== All tests passed ==='"


# Convert to Singularity if available
IMAGE_NAME="charla-nf"
VERSION="1.0.0"
DOCKER_LATEST="${IMAGE_NAME}:latest"
SINGULARITY_IMAGE="${IMAGE_NAME}_${VERSION}.sif"

singularity build "${SINGULARITY_IMAGE}" "docker-daemon://${DOCKER_LATEST}"




# Print summary
echo ""
echo "========================================"
echo "Build Summary"
echo "========================================"
echo "Docker images created:"
echo "  - ${DOCKER_TAG}"
echo "  - ${DOCKER_LATEST}"

if [ -f "${SINGULARITY_IMAGE}" ]; then
    echo ""
    echo "Singularity image created:"
    echo "  - ${SINGULARITY_IMAGE}"
    echo "  - ${IMAGE_NAME}_latest.sif (symlink)"
fi

echo ""
echo "========================================"
echo "Usage Examples"
echo "========================================"
echo ""
echo "Using Docker:"
echo "  nextflow run charla_nf/main.nf -profile docker \\"
echo "    --reads data.fq \\"
echo "    --sample_id sample1 \\"
echo "    ... other parameters ..."
echo ""
echo "Using Singularity:"
echo "  nextflow run charla_nf/main.nf -profile singularity \\"
echo "    --reads data.fq \\"
echo "    --sample_id sample1 \\"
echo "    ... other parameters ..."
echo ""

# Optional: Push to Docker Hub
echo "========================================"
echo "To push to Docker Hub:"
echo "========================================"
echo "  docker tag ${DOCKER_LATEST} yourusername/${IMAGE_NAME}:${VERSION}"
echo "  docker tag ${DOCKER_LATEST} yourusername/${IMAGE_NAME}:latest"
echo "  docker push yourusername/${IMAGE_NAME}:${VERSION}"
echo "  docker push yourusername/${IMAGE_NAME}:latest"
echo ""
