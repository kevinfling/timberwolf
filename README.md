# Timberwolf

**Time-series Impulsiveness Monitoring, Bandwise Energy Ratios, and Weighted Outlier-Limited Filtering (TIMBERWOLF).**

**Experimental condition monitoring and robust spectral feature extraction in C11.**

Timberwolf turns sampled vibration and other signals into diagnostic features and an *uncalibrated* health-index trace. It combines robust quantile statistics, spectral energy analysis, ensemble-assisted spectral decomposition through the legacy backend (with [PRISM](https://github.com/kevinfling/PRISM) designated as its replacement), and robust scalar tracking.

Originally developed as `experimental/timberwolf.h` in [CHEAP](https://github.com/kevinfling/CHEAP), Timberwolf is now a separate project. CHEAP supplies numerical operators; [PRISM](https://github.com/kevinfling/PRISM) is intended to supply decomposition; Timberwolf evaluates their utility for condition monitoring.

> **Research prototype, not a validated fault alarm.** The health index is not a calibrated failure probability, remaining-useful-life estimate, or safety-critical decision rule. Performance claims require independent experimental evidence.

## Implementation status

| Feature | Implemented behavior | Important limitation |
| --- | --- | --- |
| Robust statistics | Quantiles, MAD, Moors kurtosis, Bowley skewness | Check edge cases and median conventions |
| ASHI | DCT-II energy, Tukey weighting, cumulative-energy partitions | Adaptive percentiles may conceal spectral migration |
| Decomposition | Legacy ensemble-assisted spectral extraction; PRISM replacement pending | Not classical envelope-based CEEMDAN |
| Component weighting | Softmax of robust shape scores | Heuristic, not an optimal estimator |
| Wolf filter | Scalar robust Kalman-style update | Variance is not calibrated confidence |
| Health trace | Batch global features fused with samplewise tracker | Not independent local evidence per sample |

The current header still requires its legacy decomposition dependency and FFTW3 even when `use_spectral=0`; the PRISM integration has not yet been implemented. A dependency-free statistics/tracking core and optional PRISM/CHEAP adapters are **planned**, not shipped.

## Build and test

Requires CMake 3.18+, a C11 compiler, FFTW3 development files, Git, and network access for FetchContent.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/examples/basic_health_index
./build/benchmarks/bench_process
```

Set `-DTIMBERWOLF_BUILD_TESTS=OFF`, `-DTIMBERWOLF_BUILD_EXAMPLES=OFF`, or `-DTIMBERWOLF_BUILD_BENCHMARKS=OFF` to disable the corresponding targets. The test suite uses Unity v2.6.0. The legacy decomposition FetchContent dependency currently follows `master`; pin an exact revision for reproducible experiments. Migrating this build dependency to PRISM requires an implemented, tested API adapter.

## API

```c
#include "timberwolf.h"

size_t tw_scratch_size(size_t n);
int tw_process(const double *raw, size_t n, double temperature,
               double *health_index, double *work_buffers,
               int use_spectral);
```

`raw` and `health_index` each hold `n` doubles. Allocate **at least `tw_scratch_size(n)` bytes** for `work_buffers`: currently `4 * n * sizeof(double)`. The routine requires `n >= 8` and a finite positive temperature; it also allocates internal component storage. See [examples/basic_health_index.c](examples/basic_health_index.c) for usage. Lower-level statistics, ASHI, decomposition, and tracking routines are exposed by the header.

## Mathematical contract — current implementation

The pipeline is:

```text
signal
  -> log1p(abs(signal))
  -> global robust scale/shape statistics
  -> legacy spectral ensemble decomposition (50 members; PRISM migration pending)
  -> robust-shape softmax weighting of components
  -> weighted reconstruction
     -> global ASHI feature (optional)
     -> global robust shape feature
     -> scalar robust tracking
  -> heuristic health-index samples
```

The log transform compresses dynamic range and discards sign; it is not a general-purpose variance stabilizer.

### ASHI: adaptive spectral health index

For orthonormal DCT-II coefficients `X[k]`, the implementation robustifies energy `E[k] = X[k]^2` using Tukey weights, identifies cumulative-energy percentile boundaries, and calculates

```text
ASHI = log(1 + E_high / (E_low + epsilon)).
```

Since the boundaries are defined by *energy percentiles*, the high/low ratio may be insensitive to actual frequency movement. ASHI is an experimental feature, not a proven measure of machine degradation. Compare fixed physical-frequency bands, shaft-order bands, baseline-normalized spectra, and spectral distances.

### Ensemble decomposition

The legacy decomposition bridge maps a heuristic noise correlation parameter to a Hurst exponent and extracts spectral components. This is not a mathematically justified AR(1)-to-fGn conversion, nor does spectral sifting implement classical extrema-envelope CEEMDAN. Do not interpret these components as verified intrinsic mode functions. The planned PRISM adapter will expose explicit reconstruction contracts.

### Robust tracking and score

The Wolf filter downweights large innovations by inflating effective observation variance. It can suppress genuine abrupt faults as well as spurious outliers. The existing `1 - P` proxy is not necessarily in `[0,1]` and must not be interpreted as a probability or calibrated confidence.

The global kurtosis/skewness and ASHI terms are constant across each batch. The returned samplewise trace therefore reflects a mix of global features and the evolving filter covariance; it is **not** a truly localized anomaly measurement.

## Research and engineering roadmap

1. **Numerical audit.** Differential-test quantiles and MAD; exercise duplicates, constants, nonfinite values, short inputs, buffer boundaries, and DCT normalization. Run ASan/UBSan and strict C11 warnings.
2. **Windowed features.** Separate window-level feature extraction from temporal state tracking; add overlapping windows, timestamp/sample-rate metadata, and explicit causality.
3. **Calibrated decisions.** Fit healthy baselines and thresholds; measure false alarms per operating hour rather than assigning probabilities to heuristic scores.
4. **Decomposition ablation.** Compare direct envelope/spectral features, deterministic PRISM, and optional stochastic PRISM against a frozen historical baseline under identical budgets.
5. **Bearing-fault benchmarks.** Evaluate changing RPM/load, impulsive contamination, sensor placement, and held-out machines. Report AUROC, AUPRC, fixed-FPR sensitivity, detection delay, and full end-to-end latency.
6. **Backend modularity.** Move to reusable transform plans, caller-managed workspaces, optional CHEAP/PRISM integrations, and an independently usable C11 statistics core.

**Research hypothesis:** decomposition and robust multiscale features improve early fault detection under changing operating conditions compared with conventional envelope analysis, fixed-band spectral features, and robust time-domain statistics. This remains to be tested.

## Project layout

```text
include/timberwolf.h      experimental single-header implementation
tests/                    Unity tests
examples/                 example health-index processing
benchmarks/               end-to-end processing benchmark
CMakeLists.txt             C11 build and dependencies
```

## Reproducibility and limitations

Current tests establish some API behavior, not diagnostic validity. No detection SNR floor, speedup, calibrated alarm rate, or production readiness is claimed. Do not use `-ffast-math` for numerical correctness testing. Record compiler flags, FFT backend, dependency SHAs, seeds, signal provenance, and hardware for published results.

## License

Preserve upstream authorship and licensing obligations from the original CHEAP ecosystem. A standalone license and distribution policy should be confirmed before a tagged release.
