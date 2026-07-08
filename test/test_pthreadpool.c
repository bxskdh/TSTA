/*
 * Unit tests for pthreadpool.c -- the worker thread pool.
 *
 * These tests focus on the pool's core, deterministic behaviour: a pool
 * can be created and every task submitted with threadPoolAdd is executed
 * exactly once by a worker thread.
 *
 * IMPORTANT: threadPoolDestory neither synchronises on the pool mutex nor
 * pthread_join()s the worker threads before free()ing the pool. Calling it
 * while workers are still running (which they always are, since workers
 * block in pthread_cond_wait) is a data race / use-after-free and can hang
 * or crash nondeterministically. The tests therefore do not tear a live
 * pool down; the leaked worker threads are reaped when the test process
 * exits. Only the safe NULL-guard path of threadPoolDestory is exercised.
 * See the PR description for details of this pre-existing bug.
 *
 * The pool's worker free()s each task's arg after running it, so every task
 * argument passed to threadPoolAdd is heap allocated here.
 */

#include "minunit.h"
#include "pthreadpool.h"

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

typedef struct {
  int* counter;
  pthread_mutex_t* lock;
} taskArg;

static void
increment_task(void* arg)
{
  taskArg* a = (taskArg*)arg;
  pthread_mutex_lock(a->lock);
  (*a->counter)++;
  pthread_mutex_unlock(a->lock);
}

static int
wait_for_count(int* counter, pthread_mutex_t* lock, int target, int timeout_ms)
{
  for (int elapsed = 0; elapsed < timeout_ms; elapsed += 5) {
    pthread_mutex_lock(lock);
    int c = *counter;
    pthread_mutex_unlock(lock);
    if (c >= target) {
      return c;
    }
    usleep(5000);
  }
  pthread_mutex_lock(lock);
  int c = *counter;
  pthread_mutex_unlock(lock);
  return c;
}

static void
run_all_tasks(int worker_count, int queue_capacity, int task_count)
{
  int counter = 0;
  pthread_mutex_t lock;
  pthread_mutex_init(&lock, NULL);

  ThreadPool* pool = threadPoolCreate(worker_count, queue_capacity);
  mu_assert(pool != NULL);

  for (int i = 0; i < task_count; i++) {
    taskArg* a = (taskArg*)malloc(sizeof(taskArg));
    a->counter = &counter;
    a->lock = &lock;
    threadPoolAdd(pool, increment_task, a);
  }

  int final = wait_for_count(&counter, &lock, task_count, 5000);
  mu_assert_int_eq(task_count, final);
  /* pool intentionally not destroyed; see file header. */
}

MU_TEST(test_create_returns_pool)
{
  ThreadPool* pool = threadPoolCreate(2, 8);
  mu_assert(pool != NULL);
  /* pool intentionally not destroyed; see file header. */
}

MU_TEST(test_destroy_null_returns_error)
{
  mu_assert_int_eq(-1, threadPoolDestory(NULL));
}

MU_TEST(test_tasks_all_run_multi_worker)
{
  run_all_tasks(4, 64, 50);
}

MU_TEST(test_tasks_all_run_single_worker)
{
  run_all_tasks(1, 8, 10);
}

MU_TEST(test_more_tasks_than_queue_capacity)
{
  /* queue capacity smaller than task count exercises the producer's
   * notFull wait path in threadPoolAdd. */
  run_all_tasks(3, 4, 40);
}

int
main(void)
{
  printf("== test_pthreadpool ==\n");
  MU_RUN(test_create_returns_pool);
  MU_RUN(test_destroy_null_returns_error);
  MU_RUN(test_tasks_all_run_multi_worker);
  MU_RUN(test_tasks_all_run_single_worker);
  MU_RUN(test_more_tasks_than_queue_capacity);
  return mu_report();
}
