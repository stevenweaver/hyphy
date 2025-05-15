# HyPhy WASI Build Environment

This Docker-based build environment provides everything needed to build HyPhy with WASI threading support on an x86_64 platform.

## Setup

1. Build the Docker image:
   ```bash
   cd docker/wasi-build
   docker-compose build
   ```

2. Start the container:
   ```bash
   docker-compose up -d
   docker-compose exec wasi-build bash
   ```

3. Inside the container, run the build script:
   ```bash
   # Check out the wasi-threads branch if needed
   git checkout wasi-threads

   # Run the build script
   build-wasi.sh
   ```

## What's Included

- Ubuntu 22.04 (x86_64 architecture)
- LLVM 17 with Clang and wasm-ld
- WASI SDK 20.0
- Wasmtime (for testing)
- Build script for HyPhy WASI target with threading support

## Testing

After building, you can test the threading support:

```bash
# Run the simple threading test
wasmtime run --wasm-features=threads tests/wasi-build/threading_test.wasm

# Run the full HyPhy WASI build
wasmtime run --wasm-features=threads build/wasi-package/hyphy.wasm
```

## Artifacts

The build process creates:

1. `build/hyphy.wasm` - The WebAssembly binary file
2. `build/wasi-package/` - Directory containing the WASM file and resources
3. `tests/wasi-build/threading_test.wasm` - Simple threading test

You can copy these artifacts out of the container after building:

```bash
# From your host machine, not inside the container
docker cp <container_id>:/hyphy/build/wasi-package ./hyphy-wasi
```