#!/bin/bash
# Script to build the Docker image and start the container

set -e

cd "$(dirname "$0")"

echo "Building Docker image for WASI build environment..."
docker-compose build

echo "Starting container..."
docker-compose up -d

echo "Entering container..."
docker-compose exec wasi-build bash

# When the user exits the container shell
echo "Container session ended. You can restart with:"
echo "  cd docker/wasi-build && docker-compose exec wasi-build bash"
echo "To shut down the container:"
echo "  cd docker/wasi-build && docker-compose down"