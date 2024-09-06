#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "poa.h"
#include "../pthreadpool/pthreadpool.h"
#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>
#include <stdarg.h>
#include "seqio.h"

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


static inline void readseq(seq_s* seq, char* path) {
	seqioFastaRecord *fastaSeq = NULL;
	seqioOpenOptions opts = {
		.filename = path,
	};
	seq->seq_num = 0;
	seqioFile *file = seqioOpen(&opts);
	while ((fastaSeq = seqioReadFasta(file, fastaSeq)) != NULL) {
		seq->seq_num++;
	}
	seqioClose(file);
	fprintf(stderr, "seq_num: %d\n", seq->seq_num);
	file = seqioOpen(&opts);
	seq->S = (char**)malloc(seq->seq_num * sizeof(char*));
	int i = 0;
	while((fastaSeq = seqioReadFasta(file, fastaSeq)) != NULL) {
		size_t len = fastaSeq->sequence->length;
		seq->S[i] = (char*)malloc((len + 1) * sizeof(char));
		memcpy(seq->S[i], fastaSeq->sequence->data, len);
		seq->S[i][len] = '\0';
		i++;
	}
	fflush(stdout);
	seqioClose(file);
}

static inline void print_usage()
{
	printf("-M                      the sorce of match [default: 2]\n");
	printf("-X                      the sorce of dismatch [default: -5]\n");
	printf("-E                      the sorce of extend-gap [default: -2]\n");
	printf("-O                      the sorce of open-gap [default: -4]\n");
	printf("-T                      the number of threads [default: 10]\n");
	printf("-W                      the width of block(Multiplication of simd data width) [default: 16]\n");
	printf("-i                      the input sequence(fasta format)\n");
	printf("-o                      the output file [default: output.txt]\n");
	printf("example:\n./TSTA_msa -i seq.fa -o output.txt\n");
}

int main(int argc,char* argv[])
{
	int c;
	char* input = NULL;
	int T = 10;
	char* output = "output.txt";
	while((c = getopt(argc,argv,"M:X:E:O:T:W:i:o:")) != -1)
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
				T = atoi(optarg);
				break;
			case 'W':
				bS = atoi(optarg);
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
	L = bS * block;
	B = block;
	W = (L + B - 1) / B;
	seq_s* seq = (seq_s*)malloc(sizeof(seq_s));
	readseq(seq, input);
	int maxpthread = T;
	ThreadPool* pool = threadPoolCreate(maxpthread, 100);
	topo* s = (topo*)malloc(sizeof(topo));
	poa* p = poa_build_init(s, seq->S[0], seq->seq_num);
	s->p = p;
	for (int i = 1; i < seq->seq_num-1; i++)
	{
		s = control(s, seq->S[i], i, seq->seq_num, pool);
		s = t_sort(s, 0);
		// printf a precessBar
		if(i % 100 == 0)
		printf("\r[%d/%d]", i, seq->seq_num);
	}
	printf("\r[%d/%d]", seq->seq_num, seq->seq_num);
	printf("\n");
	s = control(s, seq->S[seq->seq_num - 1], seq->seq_num-1, seq->seq_num, pool);
	s = t_sort(s, 1);
	FILE* res = fopen(output, "w");
	printf_result(s, seq->seq_num, res);
	fclose(res);
	threadPoolDestory(pool);
	for (int i = 0; i < seq->seq_num; i++)
		free(seq->S[i]); 
	free(seq->S); seq->S = NULL;
	free(seq); seq = NULL;
	return 0;
}
