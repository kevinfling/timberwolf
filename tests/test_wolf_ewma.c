#include "unity.h"
#include "timberwolf.h"
#include <math.h>

void setUp(void) {}
void tearDown(void) {}

// ---- tw_wolf_init ----

void test_wolf_init_state(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    double scratch[8];
    tw_wolf_t wf;
    tw_wolf_init(&wf, 5.0, 1e-4, 0.1, data, 8, scratch);

    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 5.0, wf.mu);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 1.0, wf.P);
    TEST_ASSERT_EQUAL_INT(0, wf.use_hard_limit);
    TEST_ASSERT_TRUE(wf.c > 0.0);
}

void test_wolf_init_threshold_positive(void) {
    double data[16], scratch[16];
    for (int i = 0; i < 16; i++) data[i] = (double)i;
    tw_wolf_t wf;
    tw_wolf_init(&wf, 0.0, 1e-5, 0.01, data, 16, scratch);
    TEST_ASSERT_TRUE(wf.c > 0.0);
}

// ---- tw_wolf_update ----

void test_wolf_update_converges_to_constant(void) {
    // Initialize near the true value so the initial innovation is within the
    // Tukey threshold (the filter correctly rejects cold-start outliers).
    double data[32], scratch[32];
    for (int i = 0; i < 32; i++) data[i] = 3.0 + 0.01 * (i % 3 - 1);
    tw_wolf_t wf;
    tw_wolf_init(&wf, 3.0, 1e-5, 0.01, data, 32, scratch);

    // Feed 200 measurements near 3.0
    for (int i = 0; i < 200; i++) tw_wolf_update(&wf, 3.0);

    TEST_ASSERT_DOUBLE_WITHIN(0.05, 3.0, wf.mu);
}

void test_wolf_update_rejects_extreme_outlier(void) {
    double data[32], scratch[32];
    for (int i = 0; i < 32; i++) data[i] = 1.0;
    tw_wolf_t wf;
    tw_wolf_init(&wf, 1.0, 1e-5, 0.01, data, 32, scratch);

    // Converge close to 1.0 first
    for (int i = 0; i < 100; i++) tw_wolf_update(&wf, 1.0);
    double mu_before = wf.mu;

    // Inject a massive outlier
    tw_wolf_update(&wf, 1000.0);
    double mu_after = wf.mu;

    // The outlier should cause very little shift
    TEST_ASSERT_DOUBLE_WITHIN(0.1, mu_before, mu_after);
}

void test_wolf_p_decreases_with_steady_input(void) {
    double data[16], scratch[16];
    for (int i = 0; i < 16; i++) data[i] = 2.0;
    tw_wolf_t wf;
    tw_wolf_init(&wf, 2.0, 1e-6, 0.01, data, 16, scratch);

    double p0 = wf.P;
    for (int i = 0; i < 50; i++) tw_wolf_update(&wf, 2.0);

    TEST_ASSERT_TRUE(wf.P < p0);
}

int main(void) {
    UNITY_BEGIN();

    RUN_TEST(test_wolf_init_state);
    RUN_TEST(test_wolf_init_threshold_positive);
    RUN_TEST(test_wolf_update_converges_to_constant);
    RUN_TEST(test_wolf_update_rejects_extreme_outlier);
    RUN_TEST(test_wolf_p_decreases_with_steady_input);

    return UNITY_END();
}
