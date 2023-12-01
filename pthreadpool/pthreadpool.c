#include "pthreadpool.h"
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

// ����ṹ��
typedef struct Task
{
	void (*function)(void* arg);
	void* arg;
}Task;

// �̳߳ؽṹ��
struct ThreadPool
{
	// �������
	Task* taskQ;
	int queueCapacity;	//����
	int queueSize;		//��ǰ�������
	int queueFront;		//��ͷ -> ȡ����
	int queueRear;		//��β -> ������


	pthread_t managerID;	//�������߳�ID
	pthread_t* threadIDs;	//�������߳�ID
	int maxNum;				//����߳�����
	pthread_mutex_t mutexPool;	//�������̳߳�
	pthread_cond_t notFull;		//��������ǲ�������
	pthread_cond_t notEmpty;	//��������ǲ��ǿ���

	int shutdown;	//�Ƿ������̳߳أ�����Ϊ1������Ϊ0
};

ThreadPool* threadPoolCreate(int max, int queueCapacity)
{
	//�̳߳�
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

		//�������
		pool->taskQ = (Task*)malloc(sizeof(Task) * queueCapacity);
		pool->queueCapacity = queueCapacity;
		pool->queueSize = 0;
		pool->queueFront = 0;
		pool->queueRear = 0;

		pool->shutdown = 0;

		//�����߳�
		for (int i = 0; i < max; ++i)
		{
			pthread_create(&pool->threadIDs[i], NULL, worker, pool);
		}
		return pool;
	} while (0);

	//��Դ�ͷ�
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
	//�ر��̳߳�
	pool->shutdown = 1;
	//�����������������߳�
	for (int i = 0; i < pool->maxNum; ++i)
	{
		pthread_cond_signal(&pool->notEmpty);
	}
	//�ͷŶ��ڴ�
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
		//�����������߳�
		pthread_cond_wait(&pool->notFull, &pool->mutexPool);
	}
	if (pool->shutdown)
	{
		pthread_mutex_unlock(&pool->mutexPool);
		return;
	}
	//�������
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
		//��ǰ��������Ƿ�Ϊ��
		while (pool->queueSize == 0 && !pool->shutdown)
		{
			//���������߳�
			pthread_cond_wait(&pool->notEmpty, &pool->mutexPool);
		}

		//�ж��̳߳��Ƿ񱻹ر���
		if (pool->shutdown)
		{
			pthread_mutex_unlock(&pool->mutexPool);
			pthread_exit(NULL);
		}

		//�����������ȡ��һ������
		Task task;
		task.function = pool->taskQ[pool->queueFront].function;
		task.arg = pool->taskQ[pool->queueFront].arg;
		//�ƶ�ͷ���
		pool->queueFront = (pool->queueFront + 1) % pool->queueCapacity;
		pool->queueSize--;
		//����
		pthread_cond_signal(&pool->notFull);
		pthread_mutex_unlock(&pool->mutexPool);

		task.function(task.arg);
		free(task.arg);
		task.arg = NULL;

	}
	return NULL;
}
