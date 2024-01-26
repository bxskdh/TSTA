﻿#include <stdio.h>
#include <stdlib.h>
#include "pthreadpool.h"
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
int bS, L, W, M, X, E, O;
int B = block;
int* real;
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

static inline void readseq(char* argv[])
{
	FILE* fptr1 = fopen(argv[13], "r");
	FILE* fptr2 = fopen(argv[14], "r");
	int a, b;
	int offset[3];
	if (fgetc(fptr1) == '>')
	{
		fseek(fptr1, 0, 0);
		offset[0] = offset[1] = 0;
		while (fgetc(fptr1) != '\n')
			offset[0]++;
		offset[0] += 1;
		fseek(fptr1, 0, SEEK_END);
		a = ftell(fptr1) - offset[0] - 1;
		while (fgetc(fptr2) != '\n')
			offset[1]++;
		offset[1] += 1;
		fseek(fptr2, 0, SEEK_END);
		b = ftell(fptr2) - offset[1] - 1;
		if (a < b)
		{
			FILE* fptr3 = fptr2;
			fptr2 = fptr1;
			fptr1 = fptr3;
			int c = b;
			b = a;
			a = c;
			offset[2] = offset[0];
			offset[0] = offset[1];
			offset[1] = offset[2];
		}
	}
	else
	{
		fseek(fptr1, 0, SEEK_END);
		a = ftell(fptr1);
		fseek(fptr2, 0, SEEK_END);
		b = ftell(fptr2);
		if (a < b)
		{
			FILE* fptr3 = fptr2;
			fptr2 = fptr1;
			fptr1 = fptr3;
			int c = b;
			b = a;
			a = c;
		}
		offset[0] = offset[1] = 0;
	}
	fseek(fptr1, offset[0], 0);
	fseek(fptr2, offset[1], 0);
	length[0] = length[2] = a;
	length[1] = length[3] = b;
	if (a % L != 0)
		length[0] = a + (L - a % L);
	if (b % L != 0)
		length[1] = b + (L - b % L);
	//printf("1:%d,2:%d,length1:%d,length2:%d\n", a, b, length[0], length[1]);
	seq = (char**)malloc(2 * sizeof(char**));
	seq[0] = (char*)malloc((length[0] + 1) * sizeof(char));
	seq[1] = (char*)malloc((length[3] + 1) * sizeof(char));

	memset(seq[0], 'N', (size_t)length[0] + 1);
	seq[0][length[0]] = '\0';
	fgets(seq[0], a + 1, fptr1);	

	seq[1][length[3]] = '\0';
	fgets(seq[1], b + 1, fptr2);

	fclose(fptr1);
	fclose(fptr2);
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

static inline void row(int* maxsorce, int y, int pc2, char* h_s, char* t_temp, char* e_temp, char* q_temp, char* rf, int* r_temp)
{
	__mxxxi Smin, h, b2, e, f, t, s, v, mat, mis, egap, ogap, v1, h1, temp1;
	__mask mask;
	int j = 0;
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

		h1 = mm_add_epi8(s, ogap);
		f = mm_add_epi8(f, egap);
		f = mm_max_epi8(f, h1);
		f = mm_sub_epi8(f, t);

		e = mm_add_epi8(e, egap);
		e = mm_max_epi8(e, h1);
		e = mm_sub_epi8(e, v);
		mm_store(((__mxxxi*)e_temp) + x, e);

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
}

static inline void block_alignment(void* p)
{
	pack* pa = (pack*)p;
        int block_i = pa->i;
        int block_l = pa->l;

	int* maxsorce = (int*)malloc(L * sizeof(int));
	int pc0, pc2, pc4;
	char* h_s = (char*)mm_malloc(L * sizeof(char));
	char* c_temp = (char*)mm_malloc(L * sizeof(char));
        char* rf = (char*)mm_malloc(block * sizeof(char));
	
	char* t_temp = (char*)mm_malloc(L * sizeof(char));
	char* e_temp = (char*)mm_malloc(L * sizeof(char));
        char* q_temp = (char*)mm_malloc(L * sizeof(char));
	int* r_temp = (int*)mm_malloc(L * sizeof(int)); 
	
	if (block_i <= lmaxtag)
		pc0 = block_i - block_l;
	else
		pc0 = lmaxtag - block_l;
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
		row(maxsorce, i, pc2, h_s, t_temp, e_temp, q_temp, rf, r_temp);
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
	free(maxsorce);
}

int main(int argc, char* argv[])
{
        if(argc != 15)
        {
                printf("-M                      the sorce of match\n");
                printf("-X                      the sorce of dismatch\n");
                printf("-E                      the sorce of extend-gap\n");
		printf("-O                      the sorce of open-gap\n");
                printf("-T                      the number of threads\n");
                printf("-W                      the width of block(Multiplication of simd data width)\n");
                printf("[seq1] [seq2]		the input sequence(fastq)\n");
		printf("example:\n./go -M 2 -X -3 -E -2 -O -4 -T 10 -S 10 seq1.fa seq2.fa\n");
                return 0;
        }
        if(argc == 15)
        {
                M = atoi(argv[2]);
                X = atoi(argv[4]);
                E = atoi(argv[6]);
		O = atoi(argv[8]);
                bS = atoi(argv[12]);
        }
	L = block * bS;
	W = (L + B - 1) / B;

        /*if(block == 16)
                printf("support:%s\n","SSE");
        else if(block == 32)
                printf("support:%s\n","AVX");
        else
                printf("support:%s\n","AVX512");*/
	ms = MIN;

	readseq(argv);
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
	int maxpthread = atoi(argv[10]);
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
        printf("maxsorce=%d\n",ms);

	return 0;
}