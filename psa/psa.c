#include <stdio.h>
#include <stdlib.h>
#include "../pthreadpool/pthreadpool.h"
#include "seqio.h"
#include <pthread.h>
#include <string.h>
#include <unistd.h>
#include <malloc.h>
#include <stdint.h>
#if defined(__AVX512F__) && defined(__AVX512BW__) && defined(__AVX512DQ__)
#include <immintrin.h>
#elif __AVX2__
#include <immintrin.h>
#else
#include <nmmintrin.h>
#endif

#if defined(__AVX512F__) && defined(__AVX512BW__) && defined(__AVX512DQ__)
#define block 64
typedef __mmask64 __mask;
typedef __m512i __mxxxi;
#define mm_load _mm512_load_si512
#define mm_store        _mm512_store_si512
#define mm_set1_epi8    _mm512_set1_epi8
#define mm_set1_epi32    _mm512_set1_epi32
#define mm_slli(a)      _mm512_alignr_epi8(a, _mm512_inserti32x4(_mm512_shuffle_i32x4(a, a, _MM_SHUFFLE(2, 1, 0, 0)), _mm_setzero_si128(), 0), 16 - 1)
#define mm_insert_epi8(a,b,c) _mm512_inserti32x8(a, _mm256_insert_epi8(_mm512_extracti32x8_epi32(a, 0),b,c), 0)
#define mm_extract_epi8(a) _mm256_extract_epi8( _mm512_extracti32x8_epi32(a, 1),63%32)
#define mm_add_epi8     _mm512_add_epi8
#define mm_add_epi32	_mm512_add_epi32
#define mm_sub_epi8     _mm512_sub_epi8
#define mm_cmpeq_epi8   _mm512_cmpeq_epi8_mask
#define mm_max_epi8     _mm512_max_epi8
#define mm_max_epi32     _mm512_max_epi32
#define mm_setzero _mm512_setzero_si512
#define mm_blendv_epi8(a,b,c)   _mm512_mask_mov_epi8(a,c,b)
#define mm_and_epi8 _kand_mask64
#define mm_shuffle_epi32 _mm512_shuffle_epi32
#define mm_permute4x64_epi64 _mm512_shuffle_epi32
#define mm_cvtsixxx_si32 _mm512_cvtsi512_si32
//#define mm_malloc(a) aligned_alloc(64,a)
#define mm_malloc(a) aligned_malloc(a,64)
#define mm_free(a) aligned_free(a,64)
//#define mm_free(a) free(a)
#define mm0_epi8cvt32(a)       _mm512_cvtepi8_epi32(_mm512_castsi512_si128(a))
#define mm_epi8cvt32(a,b)	_mm512_cvtepi8_epi32(_mm512_extracti32x4_epi32(a, b))
#define mm_reduce_max_epi32(a) _mm512_reduce_max_epi32(a)

#elif __AVX2__
#define block 32
typedef __m256i __mask;
typedef __m256i __mxxxi;
#define mm_load _mm256_load_si256
#define mm_store        _mm256_store_si256
#define mm_set1_epi8    _mm256_set1_epi8
#define mm_set1_epi32    _mm256_set1_epi32
#define mm_slli(a)      _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1)
#define mm_insert_epi8(a,b,c)   _mm256_insert_epi8(a,b,c)
#define mm_extract_epi8(a) _mm256_extract_epi8(a,31)
#define mm_add_epi8     _mm256_add_epi8
#define mm_add_epi32	_mm256_add_epi32
#define mm_sub_epi8     _mm256_sub_epi8
#define mm_cmpeq_epi8   _mm256_cmpeq_epi8
#define mm_max_epi8     _mm256_max_epi8
#define mm_max_epi32     _mm256_max_epi32
#define mm_setzero _mm256_setzero_si256
#define mm_blendv_epi8(a,b,c)   _mm256_blendv_epi8(a,b,c)
#define mm_and_epi8 _mm256_and_si256
#define mm_shuffle_epi32 _mm256_shuffle_epi32
#define mm_permute4x64_epi64 _mm256_permute4x64_epi64
#define mm_cvtsixxx_si32 _mm256_cvtsi256_si32
//#define mm_malloc(a) aligned_alloc(32,a)
#define mm_malloc(a) aligned_malloc(a,32)
#define mm_free(a) aligned_free(a,32)
//#define mm_free(a) free(a)
#define mm0_epi8cvt32(a)	_mm256_cvtepi8_epi32(_mm256_castsi256_si128(a))
#define mm_epi8cvt32(a,b)	_mm256_cvtepi8_epi32(_mm256_castsi256_si128(_mm256_permute4x64_epi64(a, b)))
#define mm_reduce_max_epi32(a) mm256_max_reduce(a)

#else
#define block 16
typedef __m128i __mask;
typedef __m128i __mxxxi;
#define mm_load _mm_load_si128
#define mm_store        _mm_store_si128
#define mm_set1_epi8    _mm_set1_epi8
#define mm_set1_epi32    _mm_set1_epi32
#define mm_slli(a)      _mm_slli_si128(a,1)
#define mm_insert_epi8(a,b,c)   _mm_insert_epi8(a,b,c)
#define mm_extract_epi8(a) _mm_extract_epi8(a,15)
#define mm_add_epi8     _mm_add_epi8
#define mm_add_epi32	_mm_add_epi32
#define mm_sub_epi8     _mm_sub_epi8
#define mm_cmpeq_epi8   _mm_cmpeq_epi8
#define mm_max_epi8     _mm_max_epi8
#define mm_max_epi32     _mm_max_epi32
#define mm_setzero _mm_setzero_si128
#define mm_blendv_epi8(a,b,c)   _mm_blendv_epi8(a,b,c)
#define mm_and_epi8 _mm_and_si128
#define mm_shuffle_epi32 _mm_shuffle_epi32
#define mm_permute4x64_epi64 _mm_shuffle_epi32
#define mm_cvtsixxx_si32 _mm_cvtsi128_si32
#define mm_malloc(a) malloc(a)
#define mm_free(a) free(a)
#define mm0_epi8cvt32(a)	_mm_cvtepi8_epi32(a)
#define mm_epi8cvt32(a,b)	_mm_cvtepi8_epi32(_mm_srli_si128(a, b*4))
#define mm_reduce_max_epi32(a) mm128_max_reduce(a)
#endif

#define MIN -100
#define I_MIN -2000000000
#define NUM2(j) ((j) / L) * L + ((((j) % L) % W) * B + (((j) % L) / W))
int bS = 10;
int M = 2;
int X = -5;
int E = -2;
int O = -4;
int L = 10 * block;
int W = 10;
int B = block;
int* real;
char** back, ** eback, ** fback;
char* sorce, * esorce;
char* V, * F;
char** seq;
int lmaxtag, fmaxtag;
int length[4];

typedef struct pack
{
	int i;
	int l;
}pack;

volatile int ms;
pthread_mutex_t mutex;
volatile int lock;

static inline int mm128_max_reduce(__mxxxi a)
{
	__mxxxi b, c;
	b = mm_shuffle_epi32(a, _MM_SHUFFLE(3, 3, 1, 1));
	b = mm_max_epi32(a, b);
	c = mm_shuffle_epi32(b, _MM_SHUFFLE(2, 2, 2, 2));
	c = mm_max_epi32(b, c);
	return mm_cvtsixxx_si32(c);
}

static inline int mm256_max_reduce(__mxxxi e)
{
	__mxxxi f, g, h;
	f = mm_shuffle_epi32(e, _MM_SHUFFLE(3, 3, 1, 1));
	f = mm_max_epi32(e, f);
	g = mm_shuffle_epi32(f, _MM_SHUFFLE(2, 2, 2, 2));
	g = mm_max_epi32(f, g);
	h = mm_permute4x64_epi64(g, _MM_SHUFFLE(2, 2, 2, 2));
	h = mm_max_epi32(g, h);
	return mm_cvtsixxx_si32(h);
}

static inline void* aligned_malloc(size_t size, int base) {
	uint8_t* p, * q;
	//if (base <= 8) return malloc(size);
	p = malloc(size + base);
	if (p == NULL) return NULL;
	q = (uint8_t*)(((unsigned long long)(p + base)) & (~(((unsigned long long)base) - 1)));
	*(q - 1) = q - p;
	return q;
}

static inline void aligned_free(void *buffer, int base){
	uint8_t *p, *q;
	//if(base <= 8) return free(buffer);
	q = (uint8_t*)buffer;
	p = q - *(q - 1);
	free(p);
}

static inline void readseq(char* input1, char* input2)
{
	seqioFastaRecord* seq1 = NULL;
	seqioFastaRecord* seq2 = NULL;
	seqioFastaRecord* temp = NULL;
	seqioOpenOptions opts1 = {
		.filename = input1,
	};
	seqioOpenOptions opts2 = {
		.filename = input2,
	};
	seqioFile* file1 = seqioOpen(&opts1);
	seqioFile* file2 = seqioOpen(&opts2);
	seq1 = seqioReadFasta(file1, seq1);
	seq2 = seqioReadFasta(file2, seq2);
	if(seq2->sequence->length > seq1->sequence->length)
	{
		temp = seq1;
		seq1 = seq2;
		seq2 = temp;
	}
	length[0] = length[2] = seq1->sequence->length;
	length[1] = length[3] = seq2->sequence->length;
	if (seq1->sequence->length % L != 0)
		length[0] = seq1->sequence->length + (L - seq1->sequence->length % L);
	if (seq2->sequence->length % L != 0)
		length[1] = seq2->sequence->length + (L - seq2->sequence->length % L);
	seq = (char**)malloc(2 * sizeof(char**));
	seq[0] = (char*)malloc((length[0] + 1) * sizeof(char));
	seq[1] = (char*)malloc((length[3] + 1) * sizeof(char));
	memcpy(seq[0], seq1->sequence->data, seq1->sequence->length);
	seq[0][length[0]] = '\0';
	memcpy(seq[1], seq2->sequence->data, seq2->sequence->length);
	seq[1][length[3]] = '\0';
	seqioFreeRecord(seq1);
	seqioFreeRecord(seq2);
	seqioClose(file1);
	seqioClose(file2);
}

static inline void blockmatrix_init()
{
	memset(sorce, E, length[0]);
	memset(esorce, E + E + O, length[0]);
	sorce[0] = E + O;
	esorce[0] = 2 * (E + O);
	for (int i = 0; i < length[0]; i++)
		real[i] = O + ((i / L * L + ((i % L) % B) * W + ((i % L) / B)) + 1) * E;

	memset(V, E, length[3]);
	memset(F, E + E + O, length[3]);
	V[0] = E + O;
	F[0] = 2 * (E + O);
}

static inline void row(int* maxsorce, int y, int block_i, int block_l, int pc1, int pc2, int pc4, char* h_s, char* t_temp, char* e_temp, char* q_temp, char* rf, int* r_temp, char** source)
{
	__mxxxi zero, s1, s2, Smin, h, b2, e, ev, f, fv, t, s, v, mat, mis, egap, ogap, v1, h1, trace, etrace, ftrace, temp1;
	__mask mask, mask1;
	int j = 0;
	zero = mm_setzero();
	egap = mm_set1_epi8(E);
	ogap = mm_set1_epi8(O + E);
	mis = mm_set1_epi8(X);
	mat = mm_set1_epi8(M);
	Smin = mm_set1_epi8(MIN);

	f = Smin;
	f = mm_insert_epi8(f, F[pc2], 0);
	b2 = mm_set1_epi8(seq[1][pc2]);

	for (int x = 0; x < W; x++)
	{
		h = mm_load((__mxxxi*)q_temp + x);
		mask = mm_cmpeq_epi8(h, b2);
		h = mm_blendv_epi8(mis, mat, mask);
		mm_store(((__mxxxi*)h_s) + x, h);
		t = mm_load(((__mxxxi*)t_temp) + x);
		e = mm_load(((__mxxxi*)e_temp) + x);
		s = mm_max_epi8(h, e);
		s = mm_max_epi8(s, f);
		f = mm_add_epi8(f, egap);
		h1 = mm_add_epi8(s, ogap);
		f = mm_max_epi8(f, h1);
		f = mm_sub_epi8(f, t);
	}
	
	mm_store((__mxxxi*)rf, f);
	for (int x = 1; x < B; x++)
		if (rf[x - 1] + W * E - (r_temp[L - B + x] - r_temp[L - B + x - 1]) > rf[x])
			rf[x] = rf[x - 1] + W * E - (r_temp[L - B + x] - r_temp[L - B + x - 1]);
	
	f = mm_load((__mxxxi*)rf);
	temp1 = mm_sub_epi8(f, egap);
	f = mm_slli(f);
	f = mm_insert_epi8(f, F[pc2], 0);
	
	v = mm_sub_epi8(s, t);
	v = mm_max_epi8(temp1, v);
	v = mm_slli(v);
	v = mm_insert_epi8(v, V[pc2], 0);

	s1 = mm_set1_epi8(1);
	s2 = mm_add_epi8(s1, s1);
	b2 = mm_set1_epi32(I_MIN);
	for (int x = 0; x < W; x++)
	{
		h1 = mm_load(((__mxxxi*)h_s) + x);
		t = mm_load(((__mxxxi*)t_temp) + x);
		e = mm_load(((__mxxxi*)e_temp) + x);
		s = mm_max_epi8(e, f);
		s = mm_max_epi8(s, h1);
		h = mm_sub_epi8(s, v);
		mm_store(((__mxxxi*)t_temp) + x, h);

		trace = s2;
		mask = mm_cmpeq_epi8(s, f);
		trace = mm_blendv_epi8(trace, zero, mask);
		mask = mm_cmpeq_epi8(s, h1);
		trace = mm_blendv_epi8(trace, s1, mask);
		mm_store(((__mxxxi*)source[0]) + x, trace);

		h1 = mm_add_epi8(s, ogap);
		fv = mm_sub_epi8(f, v);
		f = mm_add_epi8(f, egap);
		mask1 = mm_cmpeq_epi8(f, h1);
		f = mm_max_epi8(f, h1);
		f = mm_sub_epi8(f, t);

		mask = mm_cmpeq_epi8(fv, ogap);
		ftrace = mm_blendv_epi8(s1, s2, mask);
		temp1 = mm_sub_epi8(zero, ftrace);
		mask = mm_and_epi8(mask, mask1);
		ftrace = mm_blendv_epi8(ftrace, temp1, mask);
		mm_store(((__mxxxi*)source[1]) + x, ftrace);

		ev = mm_sub_epi8(e, t);
		e = mm_add_epi8(e, egap);
		mask1 = mm_cmpeq_epi8(e, h1);
		e = mm_max_epi8(e, h1);
		e = mm_sub_epi8(e, v);
		mm_store(((__mxxxi*)e_temp) + x, e);

		mask = mm_cmpeq_epi8(ev, ogap);
		etrace = mm_blendv_epi8(s1, s2, mask);
		temp1 = mm_sub_epi8(zero, etrace);
		mask = mm_and_epi8(mask, mask1);
		etrace = mm_blendv_epi8(etrace, temp1, mask);
		mm_store(((__mxxxi*)source[2]) + x, etrace);

		v = mm_sub_epi8(s, t);
		v1 = mm0_epi8cvt32(v);
		h1 = mm_load(((__mxxxi*)r_temp) + j);
		h1 = mm_add_epi32(v1, h1);
		b2 = mm_max_epi32(b2, h1);
		mm_store(((__mxxxi*)r_temp) + j, h1);
		j++;
                v1 = mm_epi8cvt32(v, 1);
                h1 = mm_load(((__mxxxi*)r_temp) + j);
                h1 = mm_add_epi32(v1, h1);
		b2 = mm_max_epi32(b2, h1);
                mm_store(((__mxxxi*)r_temp) + j, h1); 
                j++;
                v1 = mm_epi8cvt32(v, 2);
                h1 = mm_load(((__mxxxi*)r_temp) + j);
                h1 = mm_add_epi32(v1, h1);
		b2 = mm_max_epi32(b2, h1);
                mm_store(((__mxxxi*)r_temp) + j, h1); 
                j++;
                v1 = mm_epi8cvt32(v, 3);
                h1 = mm_load(((__mxxxi*)r_temp) + j);
                h1 = mm_add_epi32(v1, h1);
		b2 = mm_max_epi32(b2, h1);
                mm_store(((__mxxxi*)r_temp) + j, h1); 
                j++;
	}
	F[pc2] = mm_extract_epi8(f);
	V[pc2] = mm_extract_epi8(v);
	maxsorce[y] = mm_reduce_max_epi32(b2);
#ifdef TRACE
	memcpy(back[pc2] + pc4, source[0], L);
        memcpy(fback[pc2] + pc4, source[1], L);
        memcpy(eback[pc2] + pc4, source[2], L);
#endif
}

static inline void block_alignment(void* p)
{
	pack* pa = (pack*)p;
        int block_i = pa->i;
        int block_l = pa->l;

	int* maxsorce = (int*)malloc(L * sizeof(int));
	int pc0, pc1, pc2, pc4;
	char* h_s = (char*)mm_malloc(L * sizeof(char));
	char* c_temp = (char*)mm_malloc(L * sizeof(char));
        char* rf = (char*)mm_malloc(block * sizeof(char));
	
	char* t_temp = (char*)mm_malloc(L * sizeof(char));
	char* e_temp = (char*)mm_malloc(L * sizeof(char));
        char* q_temp = (char*)mm_malloc(L * sizeof(char));
	int* r_temp = (int*)mm_malloc(L * sizeof(int)); 
	
	char** source = (char**)malloc(3 * sizeof(char*));
	for(int i = 0;i < 3;i++)
		source[i] = (char*)mm_malloc(L * sizeof(char));	
	if (block_i <= lmaxtag)
		pc0 = block_i - block_l;
	else
		pc0 = lmaxtag - block_l;
	pc1 = pc0 * W;
	pc4 = pc0 * L;

        memcpy(t_temp, sorce + pc4, L);
        memcpy(e_temp, esorce + pc4, L);
	memcpy(r_temp, real + pc4, L * 4);
        for(int i = 0;i < L;i++)
		q_temp[i] = seq[0][pc4 + (i % B) * W + (i / B)];
	
	for (int i = 0; i < L; i++)
	{
		if (block_i <= lmaxtag)
			pc2 = block_l * L + i;
		else
			pc2 = (block_l + block_i - lmaxtag) * L + i;
		if (pc2 >= length[3])
		{
			for(int s = i; s < L; s++)
				maxsorce[s] = I_MIN;
			break;
		}
		row(maxsorce, i, block_i, block_l, pc1, pc2, pc4, h_s, t_temp, e_temp, q_temp, rf, r_temp, source);
	}
	memcpy(sorce + pc4, t_temp, L);
        memcpy(esorce + pc4, e_temp, L);
	memcpy(real + pc4, r_temp, L * 4);

        for (int x = 1; x < L; x++)
        	if (maxsorce[0] <= maxsorce[x])
        		maxsorce[0] = maxsorce[x];

        pthread_mutex_lock(&mutex);
        if (ms <= maxsorce[0])
               ms = maxsorce[0];
        lock++;
        pthread_mutex_unlock(&mutex);
	
	mm_free(t_temp);mm_free(e_temp);mm_free(r_temp);mm_free(q_temp);
	mm_free(h_s); mm_free(rf);
	mm_free(c_temp);
	free(maxsorce);
	for(int i = 0;i < 3;i++)
		mm_free(source[i]);
	free(source);
}
#ifdef TRACE
static inline void trace(FILE* fptr)
{
	int i, j;
	int n = 0;
	i = length[3] - 1;
	j = length[2] - 1;	
	while (i >= 0 && j >= 0)
	{
		if (back[i][NUM2(j)] == 1)
		{
			i--; 
			j--;
		}
		else if (back[i][NUM2(j)] == 0)
		{
			if (j - 1 >= 0 && ((fback[i][NUM2(j)] == 1 || fback[i][NUM2(j)] == -1) || 
			((fback[i][NUM2(j)] == 2 || fback[i][NUM2(j)] == -2) && fback[i][NUM2(j - 1)] < 0)))
				back[i][NUM2(j - 1)] = 0;
			j--;
		}
		else
		{
			if (i - 1 >= 0 && ((eback[i][NUM2(j)] == 1 || eback[i][NUM2(j)] == -1) || 
			((eback[i][NUM2(j)] == 2 || eback[i][NUM2(j)] == -2) && eback[i - 1][NUM2(j)] < 0)))
				back[i - 1][NUM2(j)] = 2;
			i--;
			n++;
		}
	}
	if(i >= 0)
		n = n + i + 1;
	n = length[2] + n;
	char* a = (char*)malloc((n + 1) * sizeof(char));
	char* b = (char*)malloc((n + 1) * sizeof(char));
	a[n] = b[n] = '\0';
	i = length[3] - 1;
	j = length[2] - 1;
	while (i >= 0 && j >= 0)
	{
		if (back[i][NUM2(j)] == 1)
		{
			a[n - 1] = seq[0][j];
			b[n - 1] = seq[1][i];
			i--;
			j--;
			n--;
		}
		else if (back[i][NUM2(j)] == 0)
		{
			a[n - 1] = seq[0][j];
			b[n - 1] = '-';
			j--;
			n--;
		}
		else
		{
			a[n - 1] = '-';
			b[n - 1] = seq[1][i];
			i--;
			n--;
		}
	}
	while (j >= 0)
	{
		a[n - 1] = seq[0][j];
		b[n - 1] = '-';
		j--;
		n--;
	}
	while (i >= 0)
	{
		b[n - 1] = seq[1][i];
		a[n - 1] = '-';
		i--;
		n--;
	}

	fputs(">1\n",fptr);
	fputs(a, fptr);
	fputs("\n>2\n", fptr);
	fputs(b, fptr);
	free(a);
	free(b);
}
#endif

static inline void
print_usage(){
    printf("-M                      the sorce of match [default: 2]\n");
    printf("-X                      the sorce of dismatch [default: -5]\n");
    printf("-E                      the sorce of extend-gap [default: -2]\n");
    printf("-O                      the sorce of open-gap [default: -4]\n");
    printf("-T                      the number of threads [default: 10]\n");
    printf("-W                      the width of block(Multiplication of simd data width) [default: 16]\n");
    printf("-1                      the input sequence1(fasta)\n");
    printf("-2                      the input sequence2(fasta)\n");
#ifdef TRACE
    printf("-o                      the output file [default: output.txt]\n");
    printf("example:\n./TSTA_psa -1 seq1.fa -2 seq2.fa -o output.txt\n");
#else
	printf("example:\n./TSTA_psa -1 seq1.fa -2 seq2.fa\n");
#endif
}

int main(int argc, char* argv[])
{
    int c;
	char* input1 = NULL;
	char* input2 = NULL;
	char* output = "output.txt";
	int T = 10;
	while ((c = getopt(argc, argv, "M:X:E:O:T:W:1:2:o:")) != -1)
	{
		switch (c)
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
		case '1':
			input1 = optarg;
			break;
		case '2':
			input2 = optarg;
			break;
		case 'o':
			output = optarg;
			break;
		default:
			print_usage();
			return 0;
		}
	}

	if(input1 == NULL || input2 == NULL){
		printf("input file is not specified\n");
		print_usage();
		return 0;
	}
	L = block * bS;
	W = (L + B - 1) / B;
	ms = MIN;
	readseq(input1, input2);
	pthread_mutex_init(&mutex, NULL);
	unsigned int tsl;
	tsl = (length[0] + length[1]) / L - 1;
	fmaxtag = length[1] / L - 1;//the number of first line(the max number of block)
	lmaxtag = length[0] / L - 1;//the number of last line(the max number of block)

	sorce = (char*)malloc(length[0] * sizeof(char));//-> difference-sorce
	esorce = (char*)malloc(length[0] * sizeof(char));//↓ e-difference-sorce
	real = (int*)malloc(length[0] * sizeof(int));//real-sorce

	V = (char*)malloc(length[3] * sizeof(char));//↓ v-difference-sorce
	F = (char*)malloc(length[3] * sizeof(char));//-> f-difference-sorce
	int s;
#ifdef TRACE
	back = (char**)malloc(length[3] * sizeof(char*));
	for (int i = 0; i < length[3]; i++)
		back[i] = (char*)mm_malloc(length[0] * sizeof(char));
	eback = (char**)malloc(length[3] * sizeof(char*));
	for (int i = 0; i < length[3]; i++)
		eback[i] = (char*)mm_malloc(length[0] * sizeof(char));
	fback = (char**)malloc(length[3] * sizeof(char*));
	for (int i = 0; i < length[3]; i++)
		fback[i] = (char*)mm_malloc(length[0] * sizeof(char));
#endif
	int maxpthread = T;
	ThreadPool* pool = threadPoolCreate(maxpthread, 100);
	blockmatrix_init();
	int j = 0;
	for (int i = 0; i < tsl; i++)
	{

		if (i <= fmaxtag)
			j++;
		else if (i <= lmaxtag)
			;
		else
			j--;
		lock = 0;
		for (int l = 0; l < j; l++)
		{
			pack* pa = (pack*)malloc(sizeof(pack));
			pa->i = i;
			pa->l = l;
			threadPoolAdd(pool, block_alignment, pa);
			//block_alignment(pa);
		}
		while (lock != j) {}
	}
#ifdef TRACE
	FILE* fptr = fopen(output, "w");
	trace(fptr);
	fclose(fptr);
#endif
	pthread_mutex_destroy(&mutex);
	threadPoolDestory(pool);
	for (int i = 0; i < 2; i++)
		free(seq[i]);
	free(seq);
	free(sorce);
	free(esorce);
	free(real);
	free(F);
	free(V);
#ifdef TRACE
	for (int i = 0; i < length[3]; i++)
	{
		mm_free(back[i]);
		mm_free(eback[i]);
		mm_free(fback[i]);
	}
	free(back); 
	free(eback); 
	free(fback);
#endif
	printf("maxsorce=%d\n", ms);
	return 0;
}
