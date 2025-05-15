#ifndef _WASI_THREADING_H_
#define _WASI_THREADING_H_

/**
 * WASI threading compatibility layer for HyPhy
 * 
 * This header provides compatibility functions and definitions
 * for using pthreads with WASI, which has experimental threading
 * support.
 * 
 * The implementation supports both ARM and x86 architectures.
 */

#ifdef _USE_WASI_
#include <pthread.h>
#include <atomic>
#include <thread>

// Architecture detection (at compile time)
#if defined(__arm__) || defined(__aarch64__) || defined(_M_ARM) || defined(_M_ARM64)
  #define WASI_ARCH_ARM 1
#else
  #define WASI_ARCH_ARM 0
#endif

namespace wasi_threading {
    // Initialize WASI threading environment
    inline bool initialize() {
#if WASI_ARCH_ARM
        // ARM-specific initialization if needed
#else
        // x86 or other architecture initialization
#endif
        return true;
    }
    
    // Wrapper for pthread_create with WASI-specific handling
    inline int thread_create(pthread_t *thread, const pthread_attr_t *attr,
                           void *(*start_routine) (void *), void *arg) {
#if WASI_ARCH_ARM
        // ARM-specific thread creation handling
        // For now, we use the standard pthreads API, but this could be customized
        // if needed for ARM-specific optimizations
        return pthread_create(thread, attr, start_routine, arg);
#else
        // x86 or other architecture thread creation
        return pthread_create(thread, attr, start_routine, arg);
#endif
    }
    
    // Clean up WASI threading environment
    inline void shutdown() {
#if WASI_ARCH_ARM
        // ARM-specific cleanup
#else
        // x86 or other architecture cleanup
#endif
    }
    
    // Get optimal thread count based on hardware 
    // (this can be architecture-specific)
    inline unsigned int get_optimal_thread_count() {
        // This is a simple implementation that could be improved 
        // with more architecture-specific optimizations
        unsigned int hw_threads = std::thread::hardware_concurrency();
        
#if WASI_ARCH_ARM
        // ARM processors may have different threading characteristics
        // For now, just using a simple approach
        return hw_threads > 0 ? hw_threads : 2;
#else
        // x86 or other architectures
        return hw_threads > 0 ? hw_threads : 4;
#endif
    }
}

#endif // _USE_WASI_

#endif // _WASI_THREADING_H_