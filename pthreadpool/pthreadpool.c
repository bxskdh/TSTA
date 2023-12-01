#include "pthreadpool.h"
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

// 任务结构体
typedef struct Task
{
	void (*function)(void* arg);
	void* arg;
}Task;

// 线程池结构体
struct ThreadPool
{
	// 任务队列
	Task* taskQ;
	int queueCapacity;	//容量
	int queueSize;		//当前任务个数
	int queueFront;		//队头 -> 取数据
	int queueRear;		//队尾 -> 放数据


	pthread_t managerID;	//管理者线程ID
	pthread_t* threadIDs;	//工作的线程ID
	int maxNum;				//最大线程数量
	pthread_mutex_t mutexPool;	//锁整个线程池
	pthread_cond_t notFull;		//任务队列是不是满了
	pthread_cond_t notEmpty;	//任务队列是不是空了

	int shutdown;	//是否销毁线程池，销毁为1，否则为0
};

ThreadPool* threadPoolCreate(int max, int queueCapacity)
{
	//线程池
	ThreadPool* pool = (ThreadPool*)malloc(sizeof(ThreadPool));
	do
	{
		if (pool == NULL)
		{
			printf("malloc threadpool fail.\n");
			break;
		}
		pool->threadIDs = (pthread_t*)malloc(sizeof(pthread_t) * max);
		if (pool->threadIDs == NULL)
		{
			printf("malloc threadIDs fail.\n");
			break;
		}
		memset(pool->threadIDs, 0, sizeof(pthread_t) * max);
		pool->maxNum = max;

		if (pthread_mutex_init(&pool->mutexPool, NULL) != 0 ||
			pthread_cond_init(&pool->notEmpty, NULL) != 0 ||
			pthread_cond_init(&pool->notFull, NULL) != 0)
		{
			printf("mutex or condition init fail.\n");
			break;
		}

		//任务队列
		pool->taskQ = (Task*)malloc(sizeof(Task) * queueCapacity);
		pool->queueCapacity = queueCapacity;
		pool->queueSize = 0;
		pool->queueFront = 0;
		pool->queueRear = 0;

		pool->shutdown = 0;

		//创建线程
		for (int i = 0; i < max; ++i)
		{
			pthread_create(&pool->threadIDs[i], NULL, worker, pool);
		}
		return pool;
	} while (0);

	//资源释放
	if (pool && pool->threadIDs) free(pool->threadIDs);
	if (pool && pool->taskQ) free(pool->taskQ);
	if (pool) free(pool);

	return NULL;
}

int threadPoolDestory(ThreadPool* pool)
{
	if (pool == NULL)
	{
		return -1;
	}
	//关闭线程池
	pool->shutdown = 1;
	//唤醒阻塞的消费者线程
	for (int i = 0; i < pool->maxNum; ++i)
	{
		pthread_cond_signal(&pool->notEmpty);
	}
	//释放堆内存
	if (pool->taskQ)
	{
		free(pool->taskQ);
		pool->taskQ = NULL;
	}
	if (pool->threadIDs)
	{
		free(pool->threadIDs);
		pool->threadIDs = NULL;
	}
	pthread_mutex_destroy(&pool->mutexPool);
	pthread_cond_destroy(&pool->notEmpty);
	pthread_cond_destroy(&pool->notFull);
	free(pool);
	pool = NULL;

	return 0;
}

void threadPoolAdd(ThreadPool* pool, void(*func)(void*), void* arg)
{
	pthread_mutex_lock(&pool->mutexPool);
	while (pool->queueSize == pool->queueCapacity && !pool->shutdown)
	{
		//阻塞生产者线程
		pthread_cond_wait(&pool->notFull, &pool->mutexPool);
	}
	if (pool->shutdown)
	{
		pthread_mutex_unlock(&pool->mutexPool);
		return;
	}
	//添加任务
	pool->taskQ[pool->queueRear].function = func;
	pool->taskQ[pool->queueRear].arg = arg;
	pool->queueRear = (pool->queueRear + 1) % pool->queueCapacity;
	pool->queueSize++;

	pthread_cond_signal(&pool->notEmpty);

	pthread_mutex_unlock(&pool->mutexPool);
}

void* worker(void* arg)
{
	ThreadPool* pool = (ThreadPool*)arg;
	while (1)
	{
		pthread_mutex_lock(&pool->mutexPool);
		//当前任务队列是否为空
		while (pool->queueSize == 0 && !pool->shutdown)
		{
			//阻塞工作线程
			pthread_cond_wait(&pool->notEmpty, &pool->mutexPool);
		}

		//判断线程池是否被关闭了
		if (pool->shutdown)
		{
			pthread_mutex_unlock(&pool->mutexPool);
			pthread_exit(NULL);
		}

		//从任务队列中取出一个任务
		Task task;
		task.function = pool->taskQ[pool->queueFront].function;
		task.arg = pool->taskQ[pool->queueFront].arg;
		//移动头结点
		pool->queueFront = (pool->queueFront + 1) % pool->queueCapacity;
		pool->queueSize--;
		//解锁
		pthread_cond_signal(&pool->notFull);
		pthread_mutex_unlock(&pool->mutexPool);

		task.function(task.arg);
		free(task.arg);
		task.arg = NULL;

	}
	return NULL;
}
