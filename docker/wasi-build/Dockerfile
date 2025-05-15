FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install build essentials and dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    curl \
    unzip \
    cmake \
    python3 \
    python3-pip \
    ninja-build \
    lsb-release \
    software-properties-common \
    gnupg \
    && rm -rf /var/lib/apt/lists/*

# Install LLVM and clang
RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add - \
    && add-apt-repository "deb http://apt.llvm.org/jammy/ llvm-toolchain-jammy-17 main" \
    && apt-get update \
    && apt-get install -y \
    clang-17 \
    lld-17 \
    llvm-17 \
    && rm -rf /var/lib/apt/lists/*

# Create symbolic links for clang
RUN ln -s /usr/bin/clang-17 /usr/bin/clang \
    && ln -s /usr/bin/clang++-17 /usr/bin/clang++ \
    && ln -s /usr/bin/wasm-ld-17 /usr/bin/wasm-ld

# Download and install WASI SDK
WORKDIR /opt
RUN wget https://github.com/WebAssembly/wasi-sdk/releases/download/wasi-sdk-20/wasi-sdk-20.0-linux.tar.gz \
    && tar xf wasi-sdk-20.0-linux.tar.gz \
    && rm wasi-sdk-20.0-linux.tar.gz
ENV WASI_SDK_PATH=/opt/wasi-sdk-20.0

# Download and install Wasmtime
RUN wget https://github.com/bytecodealliance/wasmtime/releases/download/v18.0.0/wasmtime-v18.0.0-x86_64-linux.tar.xz \
    && tar xf wasmtime-v18.0.0-x86_64-linux.tar.xz \
    && mv wasmtime-v18.0.0-x86_64-linux/wasmtime /usr/local/bin/ \
    && rm -rf wasmtime-v18.0.0-x86_64-linux.tar.xz wasmtime-v18.0.0-x86_64-linux
ENV PATH="/usr/local/bin:${PATH}"

# Set up a working directory for HyPhy
WORKDIR /hyphy

# Add environment variables
ENV PATH="${WASI_SDK_PATH}/bin:${PATH}"
ENV CMAKE_TOOLCHAIN_FILE="${WASI_SDK_PATH}/share/cmake/wasi-sdk.cmake"

# Copy helper script
COPY build-wasi.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/build-wasi.sh

CMD ["/bin/bash"]