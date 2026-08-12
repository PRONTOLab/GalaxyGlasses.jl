# BayesMMfwd.jl

BayesMMfwd generates SIDES-style millimeter-source catalogs. The numerical
forward model includes star-formation rates, magnification, dust SED fluxes,
filter-grid fluxes, and CO, CII, and CI line emission.

There are two command-line programs:

| Program | Use it for | What it runs |
|---|---|---|
| `gen_pysides_from_original.jl` | Generating a catalog | Reactant only |
| `benchmark_reactant.jl` | Performance and correctness measurements | Julia and Reactant |

This separation is intentional. A normal science run should not spend time
executing the non-Reactant reference model. The benchmark runs both
implementations with identical inputs and random draws, synchronizes Reactant
before stopping its timers, and verifies every output before reporting a
speedup.

## First-time setup

All commands below are run in a terminal from the repository directory.
BayesMMfwd currently uses Julia 1.11.

Install the packages and binary artifacts once:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`--project=.` tells Julia to use this repository's `Project.toml` and
`Manifest.toml`. You should include it in every command. Each new Julia process
compiles the Reactant model for its input shape and can appear idle for roughly
one minute. The benchmark's warm samples reuse that compiled model within the
same process.

The project has been tested with Julia 1.11.7 and Reactant 0.2.259 on the CPU
backend. The local `Reactant.jl` symlink is only for Reactant developers and is
not needed for a normal science run.

### Required input files

The following files must be present at the paths specified in
`SIDES_from_original.par`:

- `SED_finegrid_dict.h5`
- `LFIR_LIR_ratio.h5`
- `Psupmu_table_Bethermin17.txt`
- `Daddi15_SLED.txt`
- a SIDES input catalog, normally under `data/`

The supplied catalogs are:

- `data/SIDES_Bethermin2017_short2.csv`: 1,024 rows for quick checks.
- `data/SIDES_Bethermin2017_short.csv`: 5,584,998 rows for production.

## Generate a catalog with Reactant

The production entrypoint is `gen_pysides_from_original.jl`. It does not run
the ordinary Julia baseline.

Start with the small catalog:

```bash
julia --project=. gen_pysides_from_original.jl
```

This compiles and executes the Reactant model, materializes its result, and
prints separate times for:

- host-to-device transfer;
- Reactant compilation;
- synchronized execution;
- device-to-host transfer;
- the total Reactant operation.

It does not write a catalog unless `--write-output` is supplied.

### Full production catalog

Before writing a catalog, open `SIDES_from_original.par` and check:

```text
output_path = "./data/SIDES/PYSIDES_ORIGINAL_OUTPUTS/"
run_name = "pySIDES_from_original"
gen_fits = true
```

Choose a new `run_name` when an existing FITS file must be preserved. Then run:

```bash
julia --project=. gen_pysides_from_original.jl --large --write-output
```

This processes all 5,584,998 rows with Reactant and writes the result using the
output settings above.

To test part of the large catalog:

```bash
julia --project=. gen_pysides_from_original.jl --large --rows 100000
```

To use another catalog:

```bash
julia --project=. gen_pysides_from_original.jl \
    --dataset path/to/catalog.csv \
    --write-output
```

The input CSV must contain `redshift`, `ra`, `dec`, `Mhalo`, and `Mstar`
columns.

Show all production options with:

```bash
julia --project=. gen_pysides_from_original.jl --help
```

### Filter-grid fluxes

Pass `--filters` to enable the filter-flux stage:

```bash
julia --project=. gen_pysides_from_original.jl \
    --large --filters --write-output
```

Each entry in `filter_list` in `SIDES_from_original.par` must have a matching
`FILTER_NAME.h5` under `grid_filter_path`. The filter grids are not currently
included in this repository. If they are absent, the driver lists every
missing file and stops. Reactant support for this stage is tested using a
temporary file-backed HDF5 grid.

## Benchmark Julia against Reactant

Use `benchmark_reactant.jl` only when collecting performance or correctness
results. It does not write a simulated catalog.

A quick correctness and timing check is:

```bash
julia --startup-file=no --history-file=no --project=. \
    benchmark_reactant.jl \
    --samples 5 \
    --results benchmark_1024.txt \
    --csv-results benchmark_1024.csv
```

For a representative 100,000-row measurement:

```bash
julia --startup-file=no --history-file=no --project=. \
    benchmark_reactant.jl \
    --large --rows 100000 \
    --samples 10 \
    --results benchmark_100k.txt \
    --csv-results benchmark_100k.csv
```

For the full catalog:

```bash
julia --startup-file=no --history-file=no --project=. \
    benchmark_reactant.jl \
    --large \
    --samples 10 \
    --results benchmark_full.txt \
    --csv-results benchmark_full.csv
```

`--startup-file=no` prevents personal Julia startup customizations from
affecting a published benchmark. `--history-file=no` is not required for
performance, but keeps noninteractive runs isolated.

Each invocation creates:

- a human-readable text report for inspection and archiving;
- a one-row CSV containing the same measurements and machine metadata for
  plotting or statistical analysis.

The report records Julia and Reactant versions, CPU model, logical CPU count,
Julia thread count, operating system, seed, catalog size, and filter count.

### Interpreting the measurements

The most useful fields are:

| Field | Meaning |
|---|---|
| `julia_warm_median_s` | Median ordinary Julia execution after Julia JIT warm-up |
| `reactant_warm_median_s` | Median synchronized execution with inputs and outputs resident in Reactant |
| `reactant_device_to_host_median_s` | Median time to copy all output arrays back to the Julia host |
| `reactant_materialized_median_s` | Median execution plus device-to-host output materialization |
| `julia_over_reactant_resident_speedup` | Julia median divided by resident Reactant execution |
| `julia_over_reactant_materialized_speedup` | Julia median divided by materialized Reactant execution |
| `reactant_compile_s` | One-time Reactant compilation cost for this input shape |
| `reactant_compile_plus_first_s` | Compilation plus the first synchronized execution |
| `reactant_one_shot_materialized_s` | Transfer, compilation, first execution, and output materialization |
| `correctness_passed` | Whether every Reactant output matched Julia |
| `maximum_relative_difference` | Largest relative difference across output arrays |

Compilation, resident execution, and materialized execution answer different
scientific-computing questions and should be reported separately. For a
one-off small catalog, compilation dominates. If later calculations remain
inside Reactant, resident execution is the appropriate throughput number. If
Julia or FITS output consumes every result, use the materialized measurement.

The numerical timers exclude CSV/HDF5 loading, luminosity-distance
preparation, DataFrame assembly, and FITS writing. Host/device transfer and
output materialization are reported explicitly.

The benchmark compares wall-clock performance on the selected host. Reactant's
CPU backend may use the CPU differently from ordinary Julia. Describe the
result as end-to-end CPU acceleration unless the processes have also been
restricted to an identical core allocation.

### Reference development results

The following runs verified the implementation on the development machine:

| Rows | Julia warm | Reactant resident | Resident speedup | Compile + first |
|---:|---:|---:|---:|---:|
| 1,024 | 0.001478 s | 0.001748 s | 0.846× | 54.406 s |
| 100,000 | 0.349141 s | 0.039502 s | 8.839× | 53.331 s |
| 5,584,998 | 27.697927 s | 1.321375 s | 20.961× | 56.013 s |

All outputs passed the Julia-versus-Reactant comparison. These numbers
demonstrate the scaling behavior, but they are machine-specific and should not
be copied directly into a paper. Rerun the benchmark protocol on the hardware
described in the manuscript.

## Recommended protocol for paper figures

For a reproducible performance figure:

1. Use a dedicated, otherwise idle machine.
2. Record the Git commit with `git rev-parse HEAD` and preserve any uncommitted
   patch used for the runs.
3. Choose catalog sizes before running the experiment.
4. Use at least 10 warm samples per size.
5. Repeat each command in at least three fresh Julia processes and retain every
   CSV file.
6. Keep the backend, hardware, thread policy, seed, filters, and input catalog
   identical across comparisons.
7. Plot Julia, Reactant resident, and Reactant materialized medians versus row
   count on logarithmic axes.
8. Show both resident and materialized speedup in a second panel or table.
9. Report Reactant compilation separately rather than adding it to every warm
   execution.
10. State that numerical equality was checked and report the tolerance or
    maximum observed difference.

On Linux, `taskset` can restrict both processes to an explicit CPU allocation.
For example, prepend `taskset -c 0-15` to the Julia command to allow CPUs 0
through 15. Use the same allocation for every run and document it.

### Example Python plotting code

The CSV files can be plotted without knowing Julia:

```python
from glob import glob

import matplotlib.pyplot as plt
import pandas as pd

results = pd.concat(
    [pd.read_csv(path) for path in glob("paper_benchmark_*.csv")],
    ignore_index=True,
).sort_values("rows")

fig, (time_ax, speedup_ax) = plt.subplots(
    2, 1, sharex=True, figsize=(6, 6),
    gridspec_kw={"height_ratios": [2, 1]},
)

time_ax.loglog(
    results["rows"], results["julia_warm_median_s"], "o-", label="Julia",
)
time_ax.loglog(
    results["rows"], results["reactant_warm_median_s"], "o-", label="Reactant",
)
time_ax.loglog(
    results["rows"],
    results["reactant_materialized_median_s"],
    "o--",
    label="Reactant + host output",
)
time_ax.set_ylabel("Warm execution time [s]")
time_ax.legend()

speedup_ax.semilogx(
    results["rows"],
    results["julia_over_reactant_resident_speedup"],
    "o-",
    label="Resident",
)
speedup_ax.semilogx(
    results["rows"],
    results["julia_over_reactant_materialized_speedup"],
    "o--",
    label="Materialized",
)
speedup_ax.axhline(1.0, color="black", linewidth=0.8)
speedup_ax.set_xlabel("Number of galaxies")
speedup_ax.set_ylabel("Speedup")
speedup_ax.legend()

fig.tight_layout()
fig.savefig("reactant_scaling.pdf")
```

For the final paper analysis, aggregate repeated-process CSV files at each row
count and add an uncertainty band or error bars rather than plotting only one
run.

## Tests

Run the regression suite after changing numerical code or input-table logic:

```bash
julia --startup-file=no --history-file=no --project=. test/runtests.jl
```

The suite compiles a complete synthetic model with Reactant, compares all
outputs with Julia, exercises filter-grid boundary cases, loads a temporary
HDF5 filter, assembles the line columns, and writes a temporary FITS file.

For a quick ordinary-Julia output smoke test:

```bash
julia --project=. gen_pysides_from_original_julia_test.jl
```

## Code layout

- `gen_pysides_from_original.jl`: Reactant-only production driver.
- `benchmark_reactant.jl`: Julia-versus-Reactant benchmark driver.
- `src/BayesMM.jl`: package entry point.
- `src/io.jl`: parameter, catalog, and FITS I/O.
- `src/reactant_pipeline.jl`: numerical model shared by Julia and Reactant.
- `src/pipeline_host.jl`: data loading, deterministic random inputs, output
  comparison, catalog assembly, timing, and reports.
- `test/runtests.jl`: regression coverage.

Numeric work intended for both backends belongs in `src/reactant_pipeline.jl`.
File I/O, DataFrames, dictionaries, and random-number generation should remain
in `src/pipeline_host.jl` so Reactant receives only arrays and simple numerical
containers. The earlier stage-specific translations were removed after their
corrected functionality and regression coverage moved into these two files.
