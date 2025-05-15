#ifndef _WASI_PARALLEL_H_
#define _WASI_PARALLEL_H_

/**
 * WASI parallel processing compatibility layer for HyPhy
 * 
 * This header provides compatibility functions for parallelizing
 * computation in WASI where OpenMP might not be available.
 * It implements simple thread pools and parallel loops to replace
 * OpenMP functionality.
 */

#include <pthread.h>
#include <vector>
#include <functional>
#include <atomic>
#include <cstddef>

#ifdef _USE_WASI_
#include <thread>

// Architecture detection (if not already defined)
#ifndef WASI_ARCH_ARM
  #if defined(__arm__) || defined(__aarch64__) || defined(_M_ARM) || defined(_M_ARM64)
    #define WASI_ARCH_ARM 1
  #else
    #define WASI_ARCH_ARM 0
  #endif
#endif

namespace wasi_parallel {

    // Simple thread pool for parallel operations
    class ThreadPool {
    public:
        ThreadPool(size_t numThreads = 0) : shutdown_(false) {
            // If numThreads is 0, determine the optimal thread count
            if (numThreads == 0) {
                // Use architecture-specific thread count
#if WASI_ARCH_ARM
                // ARM processors often benefit from fewer threads
                numThreads = std::thread::hardware_concurrency();
                numThreads = numThreads > 0 ? std::min(numThreads, 4u) : 2;
#else
                // x86 or other architectures
                numThreads = std::thread::hardware_concurrency();
                numThreads = numThreads > 0 ? numThreads : 4;
#endif
            }
            
            workers_.resize(numThreads);
            
            // Create worker threads
            for (size_t i = 0; i < numThreads; ++i) {
                pthread_create(&workers_[i], nullptr, workerThread, this);
            }
        }
        
        ~ThreadPool() {
            shutdown_ = true;
            
            // Wake up all threads to check shutdown flag
            for (size_t i = 0; i < workers_.size(); ++i) {
                pthread_cond_signal(&workCond_);
            }
            
            // Join all threads
            for (auto& thread : workers_) {
                pthread_join(thread, nullptr);
            }
        }
        
        // Execute a task for each item in the range [start, end)
        template<typename IndexType, typename Func>
        void parallel_for(IndexType start, IndexType end, Func&& func) {
            const size_t total = end - start;
            if (total <= 0) return;
            
            const size_t numThreads = workers_.size();
            if (numThreads <= 1 || total <= 1) {
                // If only one thread or a single item, just do the work directly
                for (IndexType i = start; i < end; ++i) {
                    func(i);
                }
                return;
            }
            
            // Create task package
            struct TaskInfo {
                IndexType start;
                IndexType end;
                Func* func;
                std::atomic<size_t> completed;
            };
            
            // Calculate block size for each thread
            const size_t blockSize = (total + numThreads - 1) / numThreads;
            
            // Submit tasks
            TaskInfo taskInfo = { start, end, &func, 0 };
            
            for (size_t t = 0; t < numThreads; ++t) {
                IndexType taskStart = start + t * blockSize;
                IndexType taskEnd = taskStart + blockSize;
                if (taskEnd > end) taskEnd = end;
                
                if (taskStart < taskEnd) {
                    auto task = [&taskInfo, taskStart, taskEnd]() {
                        for (IndexType i = taskStart; i < taskEnd; ++i) {
                            (*(taskInfo.func))(i);
                        }
                        taskInfo.completed.fetch_add(1);
                    };
                    
                    addTask(std::move(task));
                }
            }
            
            // Wait for all tasks to complete
            while (taskInfo.completed.load() < numThreads) {
                // Spin wait for simplicity
                // Could be replaced with condition variable
            }
        }
        
    private:
        using Task = std::function<void()>;
        std::vector<pthread_t> workers_;
        std::vector<Task> tasks_;
        pthread_mutex_t queueMutex_ = PTHREAD_MUTEX_INITIALIZER;
        pthread_cond_t workCond_ = PTHREAD_COND_INITIALIZER;
        std::atomic<bool> shutdown_;
        
        void addTask(Task&& task) {
            pthread_mutex_lock(&queueMutex_);
            tasks_.push_back(std::move(task));
            pthread_mutex_unlock(&queueMutex_);
            pthread_cond_signal(&workCond_);
        }
        
        static void* workerThread(void* arg) {
            ThreadPool* pool = static_cast<ThreadPool*>(arg);
            
            while (!pool->shutdown_) {
                Task task;
                
                {
                    pthread_mutex_lock(&pool->queueMutex_);
                    while (pool->tasks_.empty() && !pool->shutdown_) {
                        pthread_cond_wait(&pool->workCond_, &pool->queueMutex_);
                    }
                    
                    if (pool->shutdown_) {
                        pthread_mutex_unlock(&pool->queueMutex_);
                        break;
                    }
                    
                    task = std::move(pool->tasks_.back());
                    pool->tasks_.pop_back();
                    pthread_mutex_unlock(&pool->queueMutex_);
                }
                
                task();
            }
            
            return nullptr;
        }
    };
    
    // Global thread pool instance for parallel operations
    inline ThreadPool& getGlobalThreadPool() {
        static ThreadPool pool;
        return pool;
    }
    
    // Parallel for implementation to replace OpenMP parallel for
    template<typename IndexType, typename Func>
    void parallel_for(IndexType start, IndexType end, Func&& func) {
        getGlobalThreadPool().parallel_for(start, end, std::forward<Func>(func));
    }
    
} // namespace wasi_parallel
#endif // _USE_WASI_

#endif // _WASI_PARALLEL_H_