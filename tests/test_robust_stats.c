#include "unity.h"
#include "timberwolf.h"

void setUp(void) {}
void tearDown(void) {}

// ---- tw_type7_quantile ----

void test_quantile_p0_returns_first(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 1.0, tw_type7_quantile(data, 5, 0.0));
}

void test_quantile_p1_returns_last(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 5.0, tw_type7_quantile(data, 5, 1.0));
}

void test_quantile_median_odd(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 3.0, tw_type7_quantile(data, 5, 0.5));
}

void test_quantile_median_even(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0};
    // Type 7: h=(4-1)*0.5=1.5, j=1, gamma=0.5 => 0.5*2 + 0.5*3 = 2.5
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 2.5, tw_type7_quantile(data, 4, 0.5));
}

void test_quantile_single_element(void) {
    double data[] = {42.0};
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 42.0, tw_type7_quantile(data, 1, 0.5));
}

// ---- tw_robust_mad ----

void test_mad_constant_array(void) {
    double data[8], scratch[8];
    for (int i = 0; i < 8; i++) data[i] = 5.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 0.0, tw_robust_mad(data, 8, scratch));
}

void test_mad_known_sequence(void) {
    // {1,2,3,4,5}: median=3, deviations={2,1,0,1,2}, MAD=1
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double scratch[5];
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 1.0, tw_robust_mad(data, 5, scratch));
}

void test_mad_empty_returns_zero(void) {
    double scratch[1];
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 0.0, tw_robust_mad(NULL, 0, scratch));
}

// ---- tw_moors_kurtosis ----

void test_moors_kurtosis_below_8_returns_gaussian(void) {
    double data[] = {1.0, 2.0, 3.0};
    double scratch[3];
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 1.233, tw_moors_kurtosis(data, 3, scratch));
}

void test_moors_kurtosis_symmetric_near_gaussian(void) {
    // Uniform distribution — Moors kurtosis ≈ 1.168 (less than Gaussian 1.233)
    double data[64], scratch[64];
    for (int i = 0; i < 64; i++) data[i] = (double)i / 63.0;
    double k = tw_moors_kurtosis(data, 64, scratch);
    // Uniform kurtosis is less than Gaussian; check it's in a plausible range
    TEST_ASSERT_TRUE(k > 0.5 && k < 1.5);
}

// ---- tw_bowley_skew ----

void test_bowley_skew_symmetric(void) {
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double scratch[5];
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 0.0, tw_bowley_skew(data, 5, scratch));
}

void test_bowley_skew_below_4_returns_zero(void) {
    double data[] = {1.0, 2.0, 3.0};
    double scratch[3];
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 0.0, tw_bowley_skew(data, 3, scratch));
}

void test_bowley_skew_positive_skew(void) {
    // Right-skewed: more mass at the low end
    double data[] = {1.0, 1.0, 1.0, 1.0, 2.0, 5.0, 10.0, 20.0};
    double scratch[8];
    double s = tw_bowley_skew(data, 8, scratch);
    TEST_ASSERT_TRUE(s > 0.0);  // Should be positive skew
}

int main(void) {
    UNITY_BEGIN();

    RUN_TEST(test_quantile_p0_returns_first);
    RUN_TEST(test_quantile_p1_returns_last);
    RUN_TEST(test_quantile_median_odd);
    RUN_TEST(test_quantile_median_even);
    RUN_TEST(test_quantile_single_element);

    RUN_TEST(test_mad_constant_array);
    RUN_TEST(test_mad_known_sequence);
    RUN_TEST(test_mad_empty_returns_zero);

    RUN_TEST(test_moors_kurtosis_below_8_returns_gaussian);
    RUN_TEST(test_moors_kurtosis_symmetric_near_gaussian);

    RUN_TEST(test_bowley_skew_symmetric);
    RUN_TEST(test_bowley_skew_below_4_returns_zero);
    RUN_TEST(test_bowley_skew_positive_skew);

    return UNITY_END();
}
