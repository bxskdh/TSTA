#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "poa.h"
#include "../pthreadpool/pthreadpool.h"
#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>
#include <stdarg.h>

#ifdef _WIN32
#define offset 2
#elif WIN32
#define offset 2
#else
#define offset 1
#endif

char M = 2;
char X = -5;
char E = -2;
char O = -4;
int bS = 10;
int L = 10 * block;
int B = 16;
int W = 10;

typedef struct seq_s {
	int seq_num;
	char** S;
}seq_s;

seq_s* read_seq(seq_s* seq, FILE* fptr1)
{
	char last;
	int s = 0, seq_num = 0, seq_num_t = 0;
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
	free(seq_len);
	return seq;
}

static inline void print_usage()
{
	printf("-M                      the sorce of match [default: 2]\n");
	printf("-X                      the sorce of dismatch [default: -5]\n");
	printf("-E                      the sorce of extend-gap [default: -2]\n");
	printf("-O                      the sorce of open-gap [default: -4]\n");
	printf("-T                      the number of threads [default: 10]\n");
	printf("-W                      the width of block(Multiplication of simd data width) [default: 16]\n");
	printf("-i                      the input sequence(MSA/fastq)\n");
	printf("-o                      the output file [default: output.txt]\n");
	printf("example:\n./go -M 2 -X -5 -E -2 -O -4 -T 10 -S 10 -i seq.fa -f output.txt\n");
}

int main(int argc,char* argv[])
{
	int c;
	char* input = NULL;
	char* output = "output.txt";
	while((c = getopt(argc,argv,"M:X:E:O:T:W:i:f:")) != -1)
	{
		switch(c)
		{
			case 'M':
				M = atoi(optarg);
				break;
			case 'X':
				X = atoi(optarg);
				break;
			case 'E':
				E = atoi(optarg);
				break;
			case 'O':
				O = atoi(optarg);
				break;
			case 'T':
				bS = atoi(optarg);
				break;
			case 'W':
				W = atoi(optarg);
				break;
			case 'i':
				input = optarg;
				break;
			case 'o':
				output = optarg;
				break;
			default:
				print_usage();
				return 0;
		}
	}

	if(input == NULL){
		printf("input file is not specified\n");
		print_usage();
		return 0;
	}
	FILE* fptr = fopen(input, "r");
	L = bS * block;
	B = block;
	W = (L + B - 1) / B;
	seq_s* seq = (seq_s*)malloc(sizeof(seq_s));
	seq = read_seq(seq, fptr);
	fclose(fptr);
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
	FILE* res = fopen(output, "w");
	printf_result(s, seq->seq_num / 2, res);
	fclose(res);
	threadPoolDestory(pool);
	for (int i = 0; i < seq->seq_num; i++)
		free(seq->S[i]); 
	free(seq->S); seq->S = NULL;
	free(seq); seq = NULL;
	return 0;
}
