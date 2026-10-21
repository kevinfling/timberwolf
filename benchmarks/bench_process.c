#define _POSIX_C_SOURCE 199309L
#include "timberwolf.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

static double elapsed_ms(struct timespec t0, struct timespec t1) {
    return (t1.tv_sec - t0.tv_sec) * 1e3
         + (t1.tv_nsec - t0.tv_nsec) * 1e-6;
}

int main(void) {
    static const size_t sizes[] = {256, 512, 1024, 2048};
    static const int n_sizes = (int)(sizeof(sizes) / sizeof(sizes[0]));

    printf("%-10s  %-12s  %-12s\n", "n", "no-spectral", "spectral");
    printf("%-10s  %-12s  %-12s\n", "----------", "------------", "------------");

    for (int s = 0; s < n_sizes; s++) {
        size_t n = sizes[s];
        double* raw  = malloc(n * sizeof(double));
        double* hi   = malloc(n * sizeof(double));
        double* work = malloc(tw_scratch_size(n));

        for (size_t i = 0; i < n; i++)
            raw[i] = sin(2.0 * TW_PI * 5.0 * i / (double)n)
                   + 0.1 * sin(2.0 * TW_PI * 37.0 * i / (double)n);

        struct timespec t0, t1;

        clock_gettime(CLOCK_MONOTONIC, &t0);
        tw_process(raw, n, 1.0, hi, work, 0);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double ms_td = elapsed_ms(t0, t1);

        clock_gettime(CLOCK_MONOTONIC, &t0);
        tw_process(raw, n, 1.0, hi, work, 1);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double ms_sp = elapsed_ms(t0, t1);

        printf("%-10zu  %-10.3f ms  %-10.3f ms\n", n, ms_td, ms_sp);
        fflush(stdout);

        free(raw); free(hi); free(work);
    }

    return 0;
}
