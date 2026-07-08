#ifndef __minunit_h__
#define __minunit_h__

/*
 * A tiny, dependency-free unit-test harness for the TSTA C sources.
 *
 * Usage:
 *   #include "minunit.h"
 *   MU_TEST(name) { mu_assert(...); }
 *   int main(void) {
 *     MU_RUN(name);
 *     return mu_report();
 *   }
 */

#include <stdio.h>
#include <string.h>

static int mu_tests_run = 0;
static int mu_tests_failed = 0;
static int mu_checks = 0;
static int mu_current_failed = 0;

#define MU_TEST(name) static void name(void)

#define mu_check(cond, msg)                                                   \
  do {                                                                        \
    mu_checks++;                                                              \
    if (!(cond)) {                                                            \
      mu_current_failed = 1;                                                  \
      printf("    FAIL: %s:%d: %s\n", __FILE__, __LINE__, msg);               \
    }                                                                         \
  } while (0)

#define mu_assert(cond) mu_check((cond), #cond)

#define mu_assert_str_eq(expected, actual)                                    \
  do {                                                                        \
    mu_checks++;                                                              \
    if (strcmp((expected), (actual)) != 0) {                                  \
      mu_current_failed = 1;                                                  \
      printf("    FAIL: %s:%d: expected \"%s\" got \"%s\"\n", __FILE__,       \
             __LINE__, (expected), (actual));                                 \
    }                                                                         \
  } while (0)

#define mu_assert_int_eq(expected, actual)                                    \
  do {                                                                        \
    mu_checks++;                                                              \
    long _e = (long)(expected);                                               \
    long _a = (long)(actual);                                                 \
    if (_e != _a) {                                                           \
      mu_current_failed = 1;                                                  \
      printf("    FAIL: %s:%d: expected %ld got %ld\n", __FILE__, __LINE__,   \
             _e, _a);                                                         \
    }                                                                         \
  } while (0)

#define MU_RUN(test)                                                          \
  do {                                                                        \
    mu_current_failed = 0;                                                    \
    mu_tests_run++;                                                           \
    printf("  RUN  %s\n", #test);                                            \
    test();                                                                   \
    if (mu_current_failed) {                                                  \
      mu_tests_failed++;                                                      \
      printf("  ---- %s FAILED\n", #test);                                   \
    } else {                                                                  \
      printf("  ok   %s\n", #test);                                          \
    }                                                                         \
  } while (0)

static int
mu_report(void)
{
  printf("\n%d test(s), %d failed, %d checks.\n", mu_tests_run,
         mu_tests_failed, mu_checks);
  return mu_tests_failed == 0 ? 0 : 1;
}

#endif // __minunit_h__
