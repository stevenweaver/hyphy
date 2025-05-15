#include <stdio.h>
#include <pthread.h>
#include <vector>

// Simple struct to pass thread data
struct ThreadData {
    int thread_id;
    int result;
};

// Thread function
void* thread_function(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    data->result = data->thread_id * 10;
    printf("Thread %d running, result = %d\n", data->thread_id, data->result);
    return NULL;
}

int main() {
    printf("WASI Threading Test Starting...\n");
    
    const int NUM_THREADS = 4;
    std::vector<pthread_t> threads(NUM_THREADS);
    std::vector<ThreadData> thread_data(NUM_THREADS);
    
    // Create threads
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].thread_id = i;
        thread_data[i].result = 0;
        
        int ret = pthread_create(&threads[i], NULL, thread_function, &thread_data[i]);
        if (ret != 0) {
            printf("Error creating thread %d: %d\n", i, ret);
            return 1;
        }
    }
    
    // Join threads
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }
    
    // Verify results
    int sum = 0;
    for (int i = 0; i < NUM_THREADS; i++) {
        printf("Thread %d result: %d\n", i, thread_data[i].result);
        sum += thread_data[i].result;
    }
    
    printf("Sum of all results: %d\n", sum);
    printf("WASI Threading Test Complete!\n");
    
    return 0;
}