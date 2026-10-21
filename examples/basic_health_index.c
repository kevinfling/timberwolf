#include "timberwolf.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    const size_t n = 256;

    double* raw  = malloc(n * sizeof(double));
    double* hi   = malloc(n * sizeof(double));
    double* work = malloc(tw_scratch_size(n));
    if (!raw || !hi || !work) {
        fprintf(stderr, "allocation failed\n");
        return 1;
    }

    // Healthy signal: 5 Hz sine wave
    for (size_t i = 0; i < n; i++)
        raw[i] = sin(2.0 * TW_PI * 5.0 * i / (double)n);

    // Fault injection: impulse at sample 128 (simulates bearing defect)
    raw[128] += 10.0;

    fprintf(stderr, "Running tw_process (n=%zu, ensemble=%d)...\n",
            n, TW_CEEMDAN_ENSEMBLE);
    fflush(stderr);

    int rc = tw_process(raw, n, 1.0, hi, work, /*use_spectral=*/1);
    if (rc != 0) {
        fprintf(stderr, "tw_process failed\n");
        return 1;
    }
    fprintf(stderr, "Done.\n");

    printf("# sample  health_index\n");
    for (size_t i = 0; i < n; i++)
        printf("%zu\t%.6f\n", i, hi[i]);

    free(raw); free(hi); free(work);
    return 0;
}
