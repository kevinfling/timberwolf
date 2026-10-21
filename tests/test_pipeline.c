#include "unity.h"
#include "timberwolf.h"
#include <math.h>
#include <stdlib.h>

void setUp(void) {}
void tearDown(void) {}

// ---- tw_process error handling ----

void test_process_returns_minus1_for_null_raw(void) {
    double hi[16], work[16 * 4];
    TEST_ASSERT_EQUAL_INT(-1, tw_process(NULL, 16, 1.0, hi, work, 0));
}

void test_process_returns_minus1_for_null_health_index(void) {
    double raw[16], work[16 * 4];
    for (int i = 0; i < 16; i++) raw[i] = 1.0;
    TEST_ASSERT_EQUAL_INT(-1, tw_process(raw, 16, 1.0, NULL, work, 0));
}

void test_process_returns_minus1_for_null_work_buffers(void) {
    double raw[16], hi[16];
    for (int i = 0; i < 16; i++) raw[i] = 1.0;
    TEST_ASSERT_EQUAL_INT(-1, tw_process(raw, 16, 1.0, hi, NULL, 0));
}

void test_process_returns_minus1_for_n_less_than_8(void) {
    double raw[7], hi[7], work[7 * 3];
    for (int i = 0; i < 7; i++) raw[i] = 1.0;
    TEST_ASSERT_EQUAL_INT(-1, tw_process(raw, 7, 1.0, hi, work, 0));
}

// ---- tw_process output validity ----

void test_process_returns_0_for_valid_sine(void) {
    const size_t n = 256;
    double* raw  = malloc(n * sizeof(double));
    double* hi   = malloc(n * sizeof(double));
    double* work = malloc(tw_scratch_size(n));

    for (size_t i = 0; i < n; i++)
        raw[i] = sin(2.0 * TW_PI * 4.0 * i / (double)n);

    int rc = tw_process(raw, n, 1.0, hi, work, 0);
    TEST_ASSERT_EQUAL_INT(0, rc);

    free(raw); free(hi); free(work);
}

void test_process_health_index_all_finite_no_spectral(void) {
    const size_t n = 256;
    double* raw  = malloc(n * sizeof(double));
    double* hi   = malloc(n * sizeof(double));
    double* work = malloc(tw_scratch_size(n));

    for (size_t i = 0; i < n; i++)
        raw[i] = sin(2.0 * TW_PI * 4.0 * i / (double)n);

    tw_process(raw, n, 1.0, hi, work, 0);

    for (size_t i = 0; i < n; i++) {
        TEST_ASSERT_TRUE(isfinite(hi[i]));
    }

    free(raw); free(hi); free(work);
}

void test_process_health_index_all_finite_with_spectral(void) {
    const size_t n = 256;
    double* raw  = malloc(n * sizeof(double));
    double* hi   = malloc(n * sizeof(double));
    double* work = malloc(tw_scratch_size(n));

    for (size_t i = 0; i < n; i++)
        raw[i] = sin(2.0 * TW_PI * 4.0 * i / (double)n);

    tw_process(raw, n, 1.0, hi, work, 1);

    for (size_t i = 0; i < n; i++) {
        TEST_ASSERT_TRUE(isfinite(hi[i]));
    }

    free(raw); free(hi); free(work);
}

void test_process_spectral_differs_from_time_domain(void) {
    const size_t n = 256;
    double* raw   = malloc(n * sizeof(double));
    double* hi_td = malloc(n * sizeof(double));
    double* hi_sp = malloc(n * sizeof(double));
    double* work  = malloc(tw_scratch_size(n));

    for (size_t i = 0; i < n; i++)
        raw[i] = sin(2.0 * TW_PI * 4.0 * i / (double)n);

    tw_process(raw, n, 1.0, hi_td, work, 0);
    tw_process(raw, n, 1.0, hi_sp, work, 1);

    // At least one sample should differ between spectral and time-domain paths
    int any_differ = 0;
    for (size_t i = 0; i < n; i++) {
        if (fabs(hi_td[i] - hi_sp[i]) > 1e-12) { any_differ = 1; break; }
    }
    TEST_ASSERT_TRUE(any_differ);

    free(raw); free(hi_td); free(hi_sp); free(work);
}

void test_scratch_size_correct(void) {
    const size_t n = 1024;
    TEST_ASSERT_EQUAL_size_t(4 * n * sizeof(double), tw_scratch_size(n));
}

int main(void) {
    UNITY_BEGIN();

    RUN_TEST(test_process_returns_minus1_for_null_raw);
    RUN_TEST(test_process_returns_minus1_for_null_health_index);
    RUN_TEST(test_process_returns_minus1_for_null_work_buffers);
    RUN_TEST(test_process_returns_minus1_for_n_less_than_8);

    RUN_TEST(test_process_returns_0_for_valid_sine);
    RUN_TEST(test_process_health_index_all_finite_no_spectral);
    RUN_TEST(test_process_health_index_all_finite_with_spectral);
    RUN_TEST(test_process_spectral_differs_from_time_domain);
    RUN_TEST(test_scratch_size_correct);

    return UNITY_END();
}
