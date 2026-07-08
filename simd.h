#ifndef TSTA_SIMD_H
#define TSTA_SIMD_H

#include <stdint.h>
#include <stdlib.h>

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
#define mm_adds_epi8     _mm512_adds_epi8
#define mm_subs_epi8     _mm512_subs_epi8
#define mm_subs_epu8	_mm512_subs_epu8
#define mm_cmpeq_epi8   _mm512_cmpeq_epi8_mask
#define mm_cmpgt_epi8   _mm512_cmpgt_epi8_mask
#define mm_max_epi8     _mm512_max_epi8
#define mm_max_epi32     _mm512_max_epi32
#define mm_setzero _mm512_setzero_si512
#define mm_blendv_epi8(a,b,c)   _mm512_mask_mov_epi8(a,c,b)
#define mm_and_epi8 _kand_mask64
#define mm_shuffle_epi32 _mm512_shuffle_epi32
#define mm_permute4x64_epi64 _mm512_shuffle_epi32
#define mm_cvtsixxx_si32 _mm512_cvtsi512_si32
#define mm_malloc(a) aligned_malloc(a,64)
#define mm_free(a) aligned_free(a,64)
#define mm0_epi8cvt32(a)       _mm512_cvtepi8_epi32(_mm512_castsi512_si128(a))
#define mm_epi8cvt32(a,b)	_mm512_cvtepi8_epi32(_mm512_extracti32x4_epi32(a, b))
#define mm_reduce_max_epi32(a) _mm512_reduce_max_epi32(a)
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
#define mm_set1_epi32    _mm256_set1_epi32
#define mm_slli(a)      _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1)
#define mm_insert_epi8(a,b,c)   _mm256_insert_epi8(a,b,c)
#define mm_extract_epi8(a) _mm256_extract_epi8(a,31)
#define mm_add_epi8     _mm256_add_epi8
#define mm_add_epi32	_mm256_add_epi32
#define mm_sub_epi8     _mm256_sub_epi8
#define mm_adds_epi8     _mm256_adds_epi8
#define mm_subs_epi8     _mm256_subs_epi8
#define mm_subs_epu8    _mm256_subs_epu8
#define mm_cmpeq_epi8   _mm256_cmpeq_epi8
#define mm_cmpgt_epi8   _mm256_cmpgt_epi8
#define mm_max_epi8     _mm256_max_epi8
#define mm_max_epi32     _mm256_max_epi32
#define mm_setzero _mm256_setzero_si256
#define mm_blendv_epi8(a,b,c)   _mm256_blendv_epi8(a,b,c)
#define mm_and_epi8 _mm256_and_si256
#define mm_shuffle_epi32 _mm256_shuffle_epi32
#define mm_permute4x64_epi64 _mm256_permute4x64_epi64
#define mm_cvtsixxx_si32 _mm256_cvtsi256_si32
#define mm_malloc(a) aligned_malloc(a,32)
#define mm_free(a) aligned_free(a,32)
#define mm0_epi8cvt32(a)	_mm256_cvtepi8_epi32(_mm256_castsi256_si128(a))
#define mm_epi8cvt32(a,b)	_mm256_cvtepi8_epi32(_mm256_castsi256_si128(_mm256_permute4x64_epi64(a, b)))
#define mm_reduce_max_epi32(a) mm256_max_reduce(a)
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
#define mm_set1_epi32    _mm_set1_epi32
#define mm_slli(a)      _mm_slli_si128(a,1)
#define mm_insert_epi8(a,b,c)   _mm_insert_epi8(a,b,c)
#define mm_extract_epi8(a) _mm_extract_epi8(a,15)
#define mm_add_epi8     _mm_add_epi8
#define mm_add_epi32	_mm_add_epi32
#define mm_sub_epi8     _mm_sub_epi8
#define mm_adds_epi8     _mm_adds_epi8
#define mm_subs_epi8     _mm_subs_epi8
#define mm_subs_epu8    _mm_subs_epu8
#define mm_cmpeq_epi8   _mm_cmpeq_epi8
#define mm_cmpgt_epi8   _mm_cmpgt_epi8
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
#define mm_reduceto16_epi8(a) _mm_add_epi16(_mm_cvtepi8_epi16(a), _mm_cvtepi8_epi16(_mm_srli_si128(a, 8)))
#define mm_hadd_epi8(a) _mm_hadd_epi16(_mm_hadd_epi16(_mm_hadd_epi16(a, a), _mm_hadd_epi16(a, a)), _mm_hadd_epi16(_mm_hadd_epi16(a, a), _mm_hadd_epi16(a, a)))
#define mm_reduce_epi8(a) _mm_extract_epi16(a, 0)
#endif

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

#endif
