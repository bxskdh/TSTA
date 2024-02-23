#ifndef __seqio_h__
#define __seqio_h__

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define enable_gzip

#ifdef enable_gzip
#include <zlib.h>
#include <zconf.h>
#endif

#define seqioDefaultLineWidth 80
#define seqioDefaultincludeComment false
#define seqioDefaultBufferSize 1024l * 16l

#define seqioMalloc(size) malloc(size)
#define seqioRealloc(ptr, size) realloc(ptr, size)
#define seqioFree(ptr) free(ptr)

typedef enum {
  seqioRecordTypeFasta,
  seqioRecordTypeFastq,
  seqioRecordTypeUnknown
} seqioRecordType;

typedef struct {
  seqioRecordType type;
} seqioRecord;

typedef struct {
  size_t length;
  size_t capacity;
  char* data;
} seqioString;

typedef struct {
  seqioRecordType type;
  seqioString* name;
  seqioString* comment;
  seqioString* sequence;
} seqioFastaRecord;

typedef struct {
  seqioRecordType type;
  seqioString* name;
  seqioString* comment;
  seqioString* sequence;
  seqioString* quality;
} seqioFastqRecord;

typedef enum {
  seqOpenModeRead,
  seqOpenModeWrite,
  seqOpenModeAppend
} seqOpenMode;

typedef struct {
  char* filename;
  bool isGzipped;
  seqOpenMode mode;
} seqioOpenOptions;

typedef enum {
  seqioBaseCaseLower,
  seqioBaseCaseUpper,
  seqioBaseCaseOriginal
} baseCase;

typedef struct {
  size_t lineWidth;
  bool includeComment;
  baseCase baseCase;
} seqioWriteOptions;

typedef struct {
  seqioOpenOptions* options;
  void* file;
  seqioRecord* record;
  struct {
    size_t offset;
    size_t left;
    size_t capacity;
    char* data;
  } buffer;
  struct {
    seqioRecordType type;
    bool isEOF;
    int state;
  } pravite;
} seqioFile;

#define defaultSeqioWriteOptions                                              \
  {                                                                           \
    .lineWidth = seqioDefaultLineWidth,                                       \
    .includeComment = seqioDefaultincludeComment,                             \
    .baseCase = seqioBaseCaseOriginal,                                        \
  }

seqioFile* seqioOpen(seqioOpenOptions* options);
void seqioClose(seqioFile* sf);
seqioRecordType seqioGuessType(seqioFile* sf);
seqioFastaRecord* seqioReadFasta(seqioFile* sf, seqioFastaRecord* record);
seqioFastqRecord* seqioReadFastq(seqioFile* sf, seqioFastqRecord* record);
seqioRecord* seqioRead(seqioFile* sf, seqioRecord* record);
void seqioFreeRecord(void* record);
void seqioWriteFasta(seqioFile* sf,
                     seqioFastaRecord* record,
                     seqioWriteOptions* options);
void seqioWriteFastq(seqioFile* sf,
                     seqioFastqRecord* record,
                     seqioWriteOptions* options);
#endif // __seqio_h__
