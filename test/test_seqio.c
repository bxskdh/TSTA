/*
 * Unit tests for seqio.c -- the FASTA/FASTQ sequence I/O library.
 *
 * These tests exercise the public API (open/close, type guessing,
 * reading and writing of FASTA and FASTQ records, base-case conversion,
 * line wrapping and gzip round-trips), which in turn drives the internal
 * buffer, string and parsing helpers.
 */

#include "minunit.h"
#include "seqio.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static char*
make_temp_path(void)
{
  static char path[64];
  static int counter = 0;
  snprintf(path, sizeof(path), "test_seqio_tmp_%d_%d.tmp", (int)getpid(),
           counter++);
  return path;
}

static void
write_text_file(const char* path, const char* content)
{
  FILE* fp = fopen(path, "wb");
  if (fp == NULL) {
    fprintf(stderr, "cannot create temp file %s\n", path);
    exit(2);
  }
  fwrite(content, 1, strlen(content), fp);
  fclose(fp);
}

static seqioFastaRecord*
new_fasta(const char* name, const char* comment, const char* sequence)
{
  seqioFastaRecord* r =
    (seqioFastaRecord*)malloc(sizeof(seqioFastaRecord));
  r->type = seqioRecordTypeFasta;
  /* seqioStringAppend* helpers are static; build the strings via the
   * write path instead by allocating small buffers we own here. The
   * library only reads ->data/->length, so we populate them directly. */
  r->name = (seqioString*)malloc(sizeof(seqioString));
  r->comment = (seqioString*)malloc(sizeof(seqioString));
  r->sequence = (seqioString*)malloc(sizeof(seqioString));

#define FILL(field, str)                                                      \
  do {                                                                        \
    size_t _l = strlen(str);                                                  \
    (field)->capacity = _l + 1;                                              \
    (field)->length = _l;                                                    \
    (field)->data = (char*)malloc(_l + 1);                                   \
    memcpy((field)->data, (str), _l + 1);                                     \
  } while (0)

  FILL(r->name, name);
  FILL(r->comment, comment);
  FILL(r->sequence, sequence);
#undef FILL
  return r;
}

static seqioFastqRecord*
new_fastq(const char* name,
          const char* comment,
          const char* sequence,
          const char* quality)
{
  seqioFastqRecord* r =
    (seqioFastqRecord*)malloc(sizeof(seqioFastqRecord));
  r->type = seqioRecordTypeFastq;
  r->name = (seqioString*)malloc(sizeof(seqioString));
  r->comment = (seqioString*)malloc(sizeof(seqioString));
  r->sequence = (seqioString*)malloc(sizeof(seqioString));
  r->quality = (seqioString*)malloc(sizeof(seqioString));

#define FILL(field, str)                                                      \
  do {                                                                        \
    size_t _l = strlen(str);                                                  \
    (field)->capacity = _l + 1;                                              \
    (field)->length = _l;                                                    \
    (field)->data = (char*)malloc(_l + 1);                                   \
    memcpy((field)->data, (str), _l + 1);                                     \
  } while (0)

  FILL(r->name, name);
  FILL(r->comment, comment);
  FILL(r->sequence, sequence);
  FILL(r->quality, quality);
#undef FILL
  return r;
}

static void
free_manual_fasta(seqioFastaRecord* r)
{
  free(r->name->data);
  free(r->name);
  free(r->comment->data);
  free(r->comment);
  free(r->sequence->data);
  free(r->sequence);
  free(r);
}

static void
free_manual_fastq(seqioFastqRecord* r)
{
  free(r->name->data);
  free(r->name);
  free(r->comment->data);
  free(r->comment);
  free(r->sequence->data);
  free(r->sequence);
  free(r->quality->data);
  free(r->quality);
  free(r);
}

/* --------------------------------------------------------------------- */

MU_TEST(test_guess_type_fasta)
{
  char* path = make_temp_path();
  write_text_file(path, ">seq1\nACGT\n");
  seqioOpenOptions o = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* sf = seqioOpen(&o);
  mu_assert(sf != NULL);
  mu_assert_int_eq(seqioRecordTypeFasta, seqioGuessType(sf));
  seqioClose(sf);
  unlink(path);
}

MU_TEST(test_guess_type_fastq)
{
  char* path = make_temp_path();
  write_text_file(path, "@seq1\nACGT\n+\n!!!!\n");
  seqioOpenOptions o = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* sf = seqioOpen(&o);
  mu_assert(sf != NULL);
  mu_assert_int_eq(seqioRecordTypeFastq, seqioGuessType(sf));
  seqioClose(sf);
  unlink(path);
}

MU_TEST(test_read_single_fasta)
{
  char* path = make_temp_path();
  write_text_file(path, ">seq1 a description\nACGTACGT\n");
  seqioOpenOptions o = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* sf = seqioOpen(&o);
  mu_assert(sf != NULL);

  seqioFastaRecord* rec = seqioReadFasta(sf, NULL);
  mu_assert(rec != NULL);
  mu_assert_str_eq("seq1", rec->name->data);
  mu_assert_str_eq("a description", rec->comment->data);
  mu_assert_str_eq("ACGTACGT", rec->sequence->data);

  /* second read hits EOF */
  seqioFastaRecord* rec2 = seqioReadFasta(sf, rec);
  mu_assert(rec2 == NULL);

  seqioClose(sf);
  unlink(path);
}

MU_TEST(test_read_multiline_and_multi_record_fasta)
{
  char* path = make_temp_path();
  write_text_file(path,
                  ">a\nACGT\nACGT\n>b\nTTTT\n>c\nGGGGCCCC\n");
  seqioOpenOptions o = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* sf = seqioOpen(&o);
  mu_assert(sf != NULL);

  seqioFastaRecord* rec = seqioReadFasta(sf, NULL);
  mu_assert(rec != NULL);
  mu_assert_str_eq("a", rec->name->data);
  mu_assert_str_eq("ACGTACGT", rec->sequence->data);

  rec = seqioReadFasta(sf, rec);
  mu_assert(rec != NULL);
  mu_assert_str_eq("b", rec->name->data);
  mu_assert_str_eq("TTTT", rec->sequence->data);

  rec = seqioReadFasta(sf, rec);
  mu_assert(rec != NULL);
  mu_assert_str_eq("c", rec->name->data);
  mu_assert_str_eq("GGGGCCCC", rec->sequence->data);

  rec = seqioReadFasta(sf, rec);
  mu_assert(rec == NULL);

  seqioClose(sf);
  unlink(path);
}

MU_TEST(test_read_single_fastq)
{
  char* path = make_temp_path();
  write_text_file(path, "@r1 comment\nACGTA\n+\n!!!!!\n");
  seqioOpenOptions o = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* sf = seqioOpen(&o);
  mu_assert(sf != NULL);

  seqioFastqRecord* rec = seqioReadFastq(sf, NULL);
  mu_assert(rec != NULL);
  mu_assert_str_eq("r1", rec->name->data);
  mu_assert_str_eq("comment", rec->comment->data);
  mu_assert_str_eq("ACGTA", rec->sequence->data);
  mu_assert_str_eq("!!!!!", rec->quality->data);

  seqioFreeRecord(rec);
  seqioClose(sf);
  unlink(path);
}

MU_TEST(test_read_multi_record_fastq)
{
  char* path = make_temp_path();
  write_text_file(path, "@r1\nAC\n+\n##\n@r2\nGGTT\n+\nIIII\n");
  seqioOpenOptions o = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* sf = seqioOpen(&o);
  mu_assert(sf != NULL);

  seqioFastqRecord* rec = seqioReadFastq(sf, NULL);
  mu_assert(rec != NULL);
  mu_assert_str_eq("r1", rec->name->data);
  mu_assert_str_eq("AC", rec->sequence->data);
  mu_assert_str_eq("##", rec->quality->data);

  rec = seqioReadFastq(sf, rec);
  mu_assert(rec != NULL);
  mu_assert_str_eq("r2", rec->name->data);
  mu_assert_str_eq("GGTT", rec->sequence->data);
  mu_assert_str_eq("IIII", rec->quality->data);

  rec = seqioReadFastq(sf, rec);
  mu_assert(rec == NULL);

  seqioClose(sf);
  unlink(path);
}

MU_TEST(test_seqio_read_dispatch)
{
  char* path = make_temp_path();
  write_text_file(path, ">only\nACGT\n");
  seqioOpenOptions o = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* sf = seqioOpen(&o);
  mu_assert(sf != NULL);

  seqioRecord* rec = seqioRead(sf, NULL);
  mu_assert(rec != NULL);
  mu_assert_int_eq(seqioRecordTypeFasta, rec->type);
  seqioFastaRecord* fr = (seqioFastaRecord*)rec;
  mu_assert_str_eq("only", fr->name->data);
  mu_assert_str_eq("ACGT", fr->sequence->data);

  seqioFreeRecord(rec);
  seqioClose(sf);
  unlink(path);
}

MU_TEST(test_write_then_read_fasta_roundtrip)
{
  char* path = make_temp_path();
  seqioOpenOptions wo = { .filename = path, .mode = seqOpenModeWrite };
  seqioFile* wf = seqioOpen(&wo);
  mu_assert(wf != NULL);

  seqioFastaRecord* rec = new_fasta("id1", "desc", "ACGTACGTAC");
  seqioWriteOptions opt = defaultSeqioWriteOptions;
  opt.includeComment = true;
  seqioWriteFasta(wf, rec, &opt);
  seqioClose(wf);

  seqioOpenOptions ro = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* rf = seqioOpen(&ro);
  mu_assert(rf != NULL);
  mu_assert_int_eq(seqioRecordTypeFasta, seqioGuessType(rf));
  seqioFastaRecord* got = seqioReadFasta(rf, NULL);
  mu_assert(got != NULL);
  mu_assert_str_eq("id1", got->name->data);
  mu_assert_str_eq("desc", got->comment->data);
  mu_assert_str_eq("ACGTACGTAC", got->sequence->data);

  seqioFreeRecord(got);
  seqioClose(rf);
  free_manual_fasta(rec);
  unlink(path);
}

MU_TEST(test_write_fasta_linewidth_wrap)
{
  char* path = make_temp_path();
  seqioOpenOptions wo = { .filename = path, .mode = seqOpenModeWrite };
  seqioFile* wf = seqioOpen(&wo);
  mu_assert(wf != NULL);

  seqioFastaRecord* rec = new_fasta("w", "", "AAAACCCCG");
  seqioWriteOptions opt = defaultSeqioWriteOptions;
  opt.lineWidth = 4;
  seqioWriteFasta(wf, rec, &opt);
  seqioClose(wf);

  /* sequence lines should be wrapped at width 4 */
  FILE* fp = fopen(path, "rb");
  char buf[128] = { 0 };
  size_t n = fread(buf, 1, sizeof(buf) - 1, fp);
  buf[n] = '\0';
  fclose(fp);
  mu_assert_str_eq(">w\nAAAA\nCCCC\nG\n", buf);

  free_manual_fasta(rec);
  unlink(path);
}

MU_TEST(test_write_fasta_basecase_upper)
{
  char* path = make_temp_path();
  seqioOpenOptions wo = { .filename = path, .mode = seqOpenModeWrite };
  seqioFile* wf = seqioOpen(&wo);
  mu_assert(wf != NULL);

  seqioFastaRecord* rec = new_fasta("u", "", "acgt");
  seqioWriteOptions opt = defaultSeqioWriteOptions;
  opt.baseCase = seqioBaseCaseUpper;
  opt.lineWidth = 0;
  seqioWriteFasta(wf, rec, &opt);
  seqioClose(wf);

  seqioOpenOptions ro = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* rf = seqioOpen(&ro);
  seqioFastaRecord* got = seqioReadFasta(rf, NULL);
  mu_assert(got != NULL);
  mu_assert_str_eq("ACGT", got->sequence->data);
  seqioFreeRecord(got);
  seqioClose(rf);

  free_manual_fasta(rec);
  unlink(path);
}

MU_TEST(test_write_then_read_fastq_roundtrip)
{
  char* path = make_temp_path();
  seqioOpenOptions wo = { .filename = path, .mode = seqOpenModeWrite };
  seqioFile* wf = seqioOpen(&wo);
  mu_assert(wf != NULL);

  seqioFastqRecord* rec = new_fastq("q1", "cmt", "ACGTT", "IIIII");
  seqioWriteFastq(wf, rec, NULL);
  seqioClose(wf);

  seqioOpenOptions ro = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* rf = seqioOpen(&ro);
  mu_assert(rf != NULL);
  mu_assert_int_eq(seqioRecordTypeFastq, seqioGuessType(rf));
  seqioFastqRecord* got = seqioReadFastq(rf, NULL);
  mu_assert(got != NULL);
  mu_assert_str_eq("q1", got->name->data);
  mu_assert_str_eq("ACGTT", got->sequence->data);
  mu_assert_str_eq("IIIII", got->quality->data);

  seqioFreeRecord(got);
  seqioClose(rf);
  free_manual_fastq(rec);
  unlink(path);
}

MU_TEST(test_gzip_fasta_roundtrip)
{
  char* path = make_temp_path();
  seqioOpenOptions wo = {
    .filename = path, .mode = seqOpenModeWrite, .isGzipped = true
  };
  seqioFile* wf = seqioOpen(&wo);
  mu_assert(wf != NULL);
  seqioFastaRecord* rec = new_fasta("gz", "", "ACGTACGTAC");
  seqioWriteOptions opt = defaultSeqioWriteOptions;
  opt.lineWidth = 0;
  seqioWriteFasta(wf, rec, &opt);
  seqioClose(wf);

  seqioOpenOptions ro = { .filename = path, .mode = seqOpenModeRead };
  seqioFile* rf = seqioOpen(&ro);
  mu_assert(rf != NULL);
  /* magic bytes should have been detected as gzip */
  mu_assert(ro.isGzipped == true);
  mu_assert_int_eq(seqioRecordTypeFasta, seqioGuessType(rf));
  seqioFastaRecord* got = seqioReadFasta(rf, NULL);
  mu_assert(got != NULL);
  mu_assert_str_eq("gz", got->name->data);
  mu_assert_str_eq("ACGTACGTAC", got->sequence->data);

  seqioFreeRecord(got);
  seqioClose(rf);
  free_manual_fasta(rec);
  unlink(path);
}

MU_TEST(test_free_null_record_is_safe)
{
  seqioFreeRecord(NULL);
  mu_assert(1);
}

int
main(void)
{
  printf("== test_seqio ==\n");
  MU_RUN(test_guess_type_fasta);
  MU_RUN(test_guess_type_fastq);
  MU_RUN(test_read_single_fasta);
  MU_RUN(test_read_multiline_and_multi_record_fasta);
  MU_RUN(test_read_single_fastq);
  MU_RUN(test_read_multi_record_fastq);
  MU_RUN(test_seqio_read_dispatch);
  MU_RUN(test_write_then_read_fasta_roundtrip);
  MU_RUN(test_write_fasta_linewidth_wrap);
  MU_RUN(test_write_fasta_basecase_upper);
  MU_RUN(test_write_then_read_fastq_roundtrip);
  MU_RUN(test_gzip_fasta_roundtrip);
  MU_RUN(test_free_null_record_is_safe);
  return mu_report();
}
