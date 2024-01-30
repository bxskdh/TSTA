#ifndef TOPO_SORT
#define TOPO_SORT

#include <stdint.h>
#include <stdio.h>
#include "../pthreadpool/pthreadpool.h"

#define MIN -80

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
#define mm_slli(a)      _mm512_alignr_epi8(a, _mm512_inserti32x4(_mm512_shuffle_i32x4(a, a, _MM_SHUFFLE(2, 1, 0, 0)), _mm_setzero_si128(), 0), 16 - 1)
#define mm_insert_epi8(a,b,c) _mm512_inserti32x8(a, _mm256_insert_epi8(_mm512_extracti32x8_epi32(a, 0),b,c), 0)
#define mm_extract_epi8(a) _mm256_extract_epi8( _mm512_extracti32x8_epi32(a, 1),63%32)
#define mm_add_epi8     _mm512_add_epi8
#define mm_sub_epi8     _mm512_sub_epi8
#define mm_cmpeq_epi8   _mm512_cmpeq_epi8_mask
#define mm_max_epi8     _mm512_max_epi8
#define mm_setzero _mm512_setzero_si512
#define mm_blendv_epi8(a,b,c)   _mm512_mask_mov_epi8(a,c,b)
#define mm_and_epi8 _kand_mask64
#define mm_malloc(a) aligned_malloc(a,64)
#define mm_free(a) aligned_free(a,64)
#define mm_reduceto16_epi8(a) _mm512_add_epi16(_mm512_cvtepi8_epi16(_mm512_castsi512_si256(a)), _mm512_cvtepi8_epi16(_mm512_extracti32x8_epi32(a, 1)))
#define mm_hadd_epi8(a) _mm512_add_epi32(_mm512_cvtepi16_epi32(_mm512_castsi512_si256(a)), _mm512_cvtepi16_epi32(_mm512_extracti32x8_epi32(a, 1)));
#define mm_reduce_epi8(a) _mm512_reduce_add_epi32(a)

#elif __AVX2__
#define block 32
typedef __m256i __mask;
typedef __m256i __mxxxi;
#define mm_load _mm256_load_si256
#define mm_store        _mm256_store_si256
#define mm_set1_epi8    _mm256_set1_epi8
#define mm_slli(a)      _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1)
#define mm_insert_epi8(a,b,c)   _mm256_insert_epi8(a,b,c)
#define mm_extract_epi8(a) _mm256_extract_epi8(a,31)
#define mm_add_epi8     _mm256_add_epi8
#define mm_sub_epi8     _mm256_sub_epi8
#define mm_cmpeq_epi8   _mm256_cmpeq_epi8
#define mm_max_epi8     _mm256_max_epi8
#define mm_setzero _mm256_setzero_si256
#define mm_blendv_epi8(a,b,c)   _mm256_blendv_epi8(a,b,c)
#define mm_and_epi8 _mm256_and_si256
#define mm_malloc(a) aligned_malloc(a,32)
#define mm_free(a) aligned_free(a,32)
#define mm_reduceto16_epi8(a) _mm256_add_epi16(_mm256_cvtepi8_epi16(_mm256_castsi256_si128(a)), _mm256_cvtepi8_epi16(_mm256_castsi256_si128(_mm256_permute2x128_si256(a, a, 1))))
#define mm_hadd_epi8(a) _mm256_hadd_epi16(_mm256_hadd_epi16(_mm256_hadd_epi16(a, a), _mm256_hadd_epi16(a, a)), _mm256_hadd_epi16(_mm256_hadd_epi16(a, a), _mm256_hadd_epi16(a, a)))
#define mm_reduce_epi8(a) _mm256_extract_epi16(_mm256_add_epi16(a, _mm256_permute2x128_si256(a, a, 1)), 0)

#else
#define block 16
typedef __m128i __mask;
typedef __m128i __mxxxi;
#define mm_load _mm_load_si128
#define mm_store        _mm_store_si128
#define mm_set1_epi8    _mm_set1_epi8
#define mm_slli(a)      _mm_slli_si128(a,1)
#define mm_insert_epi8(a,b,c)   _mm_insert_epi8(a,b,c)
#define mm_extract_epi8(a) _mm_extract_epi8(a,15)
#define mm_add_epi8     _mm_add_epi8
#define mm_sub_epi8     _mm_sub_epi8
#define mm_cmpeq_epi8   _mm_cmpeq_epi8
#define mm_max_epi8     _mm_max_epi8
#define mm_setzero _mm_setzero_si128
#define mm_blendv_epi8(a,b,c)   _mm_blendv_epi8(a,b,c)
#define mm_and_epi8 _mm_and_si128
#define mm_malloc(a) malloc(a)
#define mm_free(a) free(a)
#define mm_reduceto16_epi8(a) _mm_add_epi16(_mm_cvtepi8_epi16(a), _mm_cvtepi8_epi16(_mm_srli_si128(a, 8)))
#define mm_hadd_epi8(a) _mm_hadd_epi16(_mm_hadd_epi16(_mm_hadd_epi16(a, a), _mm_hadd_epi16(a, a)), _mm_hadd_epi16(_mm_hadd_epi16(a, a), _mm_hadd_epi16(a, a)))
#define mm_reduce_epi8(a) _mm_extract_epi16(a, 0)
#endif

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

topo* t_sort(topo* g, int num);

//poa
poa* poa_build_init(topo* p, char a[], int sum);
topo* control(topo* p, char* A, int num, int sum, ThreadPool* pool);
void printf_result(topo* p, int num, FILE* res);
#endif
