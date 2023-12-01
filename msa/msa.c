#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "poa.h"
#include "../pthreadpool/pthreadpool.h"
#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>

#ifdef _WIN32
#define offset 2
#elif WIN32
#define offset 2
#else
#define offset 1
#endif

char M, X, E, O;
int bS, L, B, W;

typedef struct seq_s {
	int seq_num;
	char** S;
}seq_s;

seq_s* read_seq(seq_s* seq, char* argv[])
{
	char last;
	int s = 0, seq_num = 0, seq_num_t = 0;
	FILE* fptr1 = fopen(argv[14], "r");
	while (!feof(fptr1))
		if (fgetc(fptr1) == '\n')
			seq_num++;
	fseek(fptr1, -1, 2);
	last = fgetc(fptr1);
	if (last != '\n')
		seq_num = seq_num + 1; 
	seq->seq_num = seq_num;
	int* seq_len = (int*)malloc(seq_num * sizeof(int));
	fseek(fptr1, 0, 0);
	while (!feof(fptr1))
	{
		if (fgetc(fptr1) == '\n')
		{
			seq_len[seq_num_t] = s;
			s = 0;
			seq_num_t++;
			continue;
		}
		s++;
	}
	if (s != 1)
		seq_len[seq_num - 1] = s - 1;
	seq->S = (char**)malloc(seq_num * sizeof(char*));
	for (int i = 0; i < seq_num; i++)
		seq->S[i] = (char*)malloc((seq_len[i] + 1) * sizeof(char));
	fseek(fptr1, 0, 0);
	for (int i = 0; i < seq_num; i++)
	{
		seq->S[i][seq_len[i]] = '\0';
		fgets(seq->S[i], seq_len[i] + 1, fptr1);
		fseek(fptr1, offset, 1);
		//printf("%d\n", strlen(seq->S[i]));
	}
	fclose(fptr1);
	free(seq_len);
	return seq;
}

int main(int argc,char* argv[])
{
	if(argc != 15)
        {
                printf("-M                      the sorce of match\n");
                printf("-X                      the sorce of dismatch\n");
                printf("-E                      the sorce of extend-gap\n");
		printf("-O                      the sorce of open-gap\n");
                printf("-T                      the number of threads\n");
                printf("-W                      the width of block(Multiplication of simd data width)\n");
                printf("-i                      the input sequence(MSA/fastq)\n");
                printf("example:\n./go -M 2 -X -5 -E -2 -O -4 -T 10 -S 10 -i seq.fa\n");
                return 0;
        }
        else
        {
                M = atoi(argv[2]);
                X = atoi(argv[4]);
                E = atoi(argv[6]);
		O = atoi(argv[8]);
                bS = atoi(argv[12]);
        }
        if(block == 16)
		printf("support:%s\n","SSE");
        else if(block == 32)
		printf("support:%s\n","AVX");
        else
		printf("support:%s\n","AVX512");
	L = bS * block;
	B = block;
	W = (L + B - 1) / B;
	
	seq_s* seq = (seq_s*)malloc(sizeof(seq_s));
	seq = read_seq(seq, argv);
	int maxpthread = atoi(argv[10]);
	ThreadPool* pool = threadPoolCreate(maxpthread, 100);
	topo* s = (topo*)malloc(sizeof(topo));
	poa* p = poa_build_init(s, seq->S[1], seq->seq_num / 2);
	s->p = p;
	for (int i = 1; i < seq->seq_num / 2 - 1; i++)
	{
		s = control(s, seq->S[i * 2 + 1], i, seq->seq_num / 2, pool);
		s = t_sort(s, 0);
	}
	s = control(s, seq->S[seq->seq_num - 1], seq->seq_num / 2 - 1, seq->seq_num / 2, pool);
	s = t_sort(s, 1);
	printf_result(s, seq->seq_num / 2);
	threadPoolDestory(pool);	//´Ý»ÙÏß³Ì³Ø
	for (int i = 0; i < seq->seq_num; i++)
		free(seq->S[i]); 
	free(seq->S); seq->S = NULL;
	free(seq); seq = NULL;
	return 0;
}
