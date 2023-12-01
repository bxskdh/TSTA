#ifndef _THREADPOOL_H
#define _THREADPOOL_H

typedef struct ThreadPool ThreadPool;
// �����̳߳ز���ʼ��
ThreadPool* threadPoolCreate(int max, int queueCapacity);
//ThreadPool* threadPoolCreate(int min, int max, int queueCapacity);

// �����̳߳�
int threadPoolDestory(ThreadPool* pool);

// ���̳߳��������
void threadPoolAdd(ThreadPool* pool, void(*func)(void*), void* arg);

/////////////////////
void* worker(void* arg);

// ��ȡ�̳߳��й������̸߳���
//int threadPoolBusyNum(ThreadPool* pool);
//int threadPoolAliveNum(ThreadPool* pool);
//void* worker(void* arg);
//void* manager(void* arg);
//void threadExit(ThreadPool* pool);

#endif // _THREADPOOL_H
