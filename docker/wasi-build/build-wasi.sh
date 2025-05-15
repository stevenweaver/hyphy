#!/bin/bash
# Script to build HyPhy WASI target with threading support

set -e

# Working directory is /hyphy

# Clean build
rm -rf build
mkdir -p build
cd build

# Configure for WASI
cmake -DCMAKE_TOOLCHAIN_FILE="${WASI_SDK_PATH}/share/cmake/wasi-sdk.cmake" \
      -DWASI_SDK_PREFIX="${WASI_SDK_PATH}" \
      -DCMAKE_C_COMPILER="${WASI_SDK_PATH}/bin/clang" \
      -DCMAKE_CXX_COMPILER="${WASI_SDK_PATH}/bin/clang++" \
      -DCMAKE_C_FLAGS="-pthread" \
      -DCMAKE_CXX_FLAGS="-pthread" \
      -DCMAKE_EXE_LINKER_FLAGS="-pthread -Wl,--import-memory,--export-memory,--shared-memory,--max-memory=1073741824 -Wl,--export-all" \
      -DCMAKE_SYSTEM_NAME=WASI \
      -DCMAKE_SYSTEM_VERSION=1 \
      -DCMAKE_SYSTEM_PROCESSOR=wasm32 \
      -DCMAKE_SYSROOT="${WASI_SDK_PATH}/share/wasi-sysroot" \
      ..

# Build
make -j$(nproc) hyphy

echo "Build completed! Artifacts:"
ls -la

# Create a package directory
mkdir -p wasi-package
cp hyphy.wasm wasi-package/
cp -r ../res wasi-package/

echo "WASI package created at build/wasi-package"
echo "To run with Wasmtime, use: wasmtime run --wasm-features=threads wasi-package/hyphy.wasm"

# Also build and run a simple threading test
cd ..
mkdir -p tests/wasi-build
cat > tests/wasi-build/threading_test.c << EOL
#include <stdio.h>
#include <pthread.h>

#define NUM_THREADS 4
#define ITERATIONS 1000

// Shared counter (using mutex for thread safety)
int counter = 0;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

// Thread function that increments the counter
void* increment_counter(void* arg) {
    int thread_id = *((int*)arg);
    
    for (int i = 0; i < ITERATIONS; i++) {
        pthread_mutex_lock(&mutex);
        counter++;
        pthread_mutex_unlock(&mutex);
    }
    
    printf("Thread %d completed\\n", thread_id);
    return NULL;
}

int main() {
    printf("WASI Threading Test Starting...\\n");
    
    pthread_t threads[NUM_THREADS];
    int thread_ids[NUM_THREADS];
    
    // Create threads
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_ids[i] = i;
        int ret = pthread_create(&threads[i], NULL, increment_counter, &thread_ids[i]);
        if (ret != 0) {
            printf("Error creating thread %d: %d\\n", i, ret);
            return 1;
        }
        printf("Created thread %d\\n", i);
    }
    
    // Join threads
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
        printf("Joined thread %d\\n", i);
    }
    
    // Verify final counter value
    printf("Final counter value: %d (expected: %d)\\n", counter, NUM_THREADS * ITERATIONS);
    printf("WASI Threading Test Complete!\\n");
    
    return 0;
}
EOL

# Build the test
"${WASI_SDK_PATH}/bin/clang" \
    --target=wasm32-wasi-threads \
    --sysroot="${WASI_SDK_PATH}/share/wasi-sysroot" \
    -pthread \
    -Wl,--import-memory,--export-memory,--shared-memory,--max-memory=67108864 \
    -O2 \
    -o tests/wasi-build/threading_test.wasm \
    tests/wasi-build/threading_test.c

echo "Threading test built at tests/wasi-build/threading_test.wasm"
echo "To run the test: wasmtime run --wasm-features=threads tests/wasi-build/threading_test.wasm"