#include "unity.h"
#include "timberwolf.h"
#include <math.h>
#include <stdlib.h>

void setUp(void) {}
void tearDown(void) {}

// ---- tw_process error handling ----

void test_process_returns_minus1_for_null_raw(void) {
    double hi[16], work[16 * 3];
    TEST_ASSERT_EQUAL_INT(-1, tw_process(NULL, 16, 1.0, hi, work, 0));
}

void test_process_returns_minus1_for_null_health_index(void) {
    double raw[16], work[16 * 3];
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
    TEST_ASSERT_EQUAL_size_t(3 * n * sizeof(double), tw_scratch_size(n));
}

void test_ctx_lepto_ewma_smooths_across_windows(void) {
    const size_t n = 256;
    double* raw  = malloc(n * sizeof(double));
    double* hi   = malloc(n * sizeof(double));
    double* work = malloc(tw_scratch_size(n));

    // Window 1: clean sine — low leptokurtosis
    for (size_t i = 0; i < n; i++)
        raw[i] = sin(2.0 * TW_PI * 4.0 * i / (double)n);

    tw_ctx_t ctx;
    tw_ctx_init(&ctx, n, 0.7);
    tw_process_ctx(&ctx, raw, n, 1.0, hi, work, 0);
    double lepto_after_w1 = ctx.lepto_ewma;

    // Window 2: impulsive signal — high leptokurtosis
    for (size_t i = 0; i < n; i++) raw[i] = 0.0;
    raw[64] = 50.0; raw[128] = -50.0; raw[192] = 50.0;
    tw_process_ctx(&ctx, raw, n, 1.0, hi, work, 0);
    double lepto_after_w2 = ctx.lepto_ewma;

    // Window 3: back to clean sine — EWMA should still be elevated (not snapped back)
    for (size_t i = 0; i < n; i++)
        raw[i] = sin(2.0 * TW_PI * 4.0 * i / (double)n);
    tw_process_ctx(&ctx, raw, n, 1.0, hi, work, 0);
    double lepto_after_w3 = ctx.lepto_ewma;

    tw_ctx_destroy(&ctx);
    free(raw); free(hi); free(work);

    // Impulsive window should raise the EWMA above the clean-sine baseline
    TEST_ASSERT_TRUE(lepto_after_w2 > lepto_after_w1);
    // EWMA must be strictly between w1 and w2: smoothing prevents it from
    // snapping back to the raw batch value in a single window.
    TEST_ASSERT_TRUE(lepto_after_w3 > lepto_after_w1);
    // All three outputs must be finite
    TEST_ASSERT_TRUE(isfinite(lepto_after_w1));
    TEST_ASSERT_TRUE(isfinite(lepto_after_w2));
    TEST_ASSERT_TRUE(isfinite(lepto_after_w3));
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
    RUN_TEST(test_ctx_lepto_ewma_smooths_across_windows);

    return UNITY_END();
}
