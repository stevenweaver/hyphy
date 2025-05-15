#!/bin/bash
# Build script for WASI threading test

# Check for WASI_SDK_PATH environment variable
if [ -z "$WASI_SDK_PATH" ]; then
    echo "WASI_SDK_PATH is not set. Please set it to the path of your WASI SDK installation."
    exit 1
fi

echo "Building threading test with WASI SDK at $WASI_SDK_PATH"

# Create build directory
mkdir -p build

# Compile the test
$WASI_SDK_PATH/bin/clang++ \
    --target=wasm32-wasi-threads \
    -pthread \
    -Wl,--import-memory,--export-memory,--shared-memory,--max-memory=67108864 \
    -O2 \
    -o build/threading_test.wasm \
    threading_test.cpp

echo "Build complete. Test executable is at build/threading_test.wasm"
echo "To run with Wasmtime:"
echo "wasmtime run --wasm-features=threads build/threading_test.wasm"