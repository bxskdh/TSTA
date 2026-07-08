#ifndef TOPO_SORT
#define TOPO_SORT

#include <stdint.h>
#include <stdio.h>
#include "../pthreadpool/pthreadpool.h"
#include "../simd.h"

#define MIN -120

typedef struct poa {
	struct poa** pre;
	struct poa** next;
	char* sorce;
	char* esorce;
	char* source;
	char* esource;
	char* fsource;
	char* passing_seq;
	struct poa* mismatch_node[4];
	int mismatch_num;
	int* simple_sorce;
	int sub;
	int frist_col_sorce;
	int in_temp;
	int in;
	int out;
	char base;
	char* f0;
	int lastsorce;
	int node_logo;
	int passing;
	int node_sorce;
	int node_sorce_source;
	int node_base_len;
	int* edge_weight;
}poa;

typedef struct topo {
	int len;
	int last_node_num;
	poa* p;
	poa** unsort;
	poa** sort;
}topo;

topo* t_sort(topo* g, int num);

//poa
poa* poa_build_init(topo* p, char a[], int sum);
topo* control(topo* p, char* A, int num, int sum, ThreadPool* pool);
void printf_result(topo* p, int num, FILE* res);
#endif
