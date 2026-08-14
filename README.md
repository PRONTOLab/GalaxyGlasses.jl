# BayesMMfwd.jl

BayesMMfwd.jl implements the numerical forward model used to generate
SIDES-style millimeter-source catalogs. Given galaxy properties, model
parameters, and explicit random draws, it computes:

- quenched and starburst classifications and star-formation rates;
- gravitational magnification;
- dust SED quantities, monochromatic fluxes, and optional filter fluxes;
- CO, [CII], and [CI] line luminosities and fluxes.

The core entry point is `BayesMMfwd.forward_model(inputs, params, noise)`. It is
pure numerical Julia: file loading, random-number generation, backend
selection, compilation, and catalog output remain outside the function. The
same forward-model source therefore runs as ordinary Julia or compiles with
Reactant for CPU, NVIDIA GPU, and TPU execution.

## Setup

BayesMMfwd currently targets Julia 1.11. From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The paths in `data/SIDES_from_original.par` must resolve to the model tables:

- `SED_finegrid_dict.h5`
- `LFIR_LIR_ratio.h5`
- `data/Psupmu_table_Bethermin17.txt`
- `data/Daddi15_SLED.txt`
- a SIDES CSV catalog

Input catalogs must contain `redshift`, `Mstar`, `Mhalo`, `ra`, and
`dec`. Filter-grid HDF5 files are required only when filter fluxes are
enabled.

## Forward-model API

```julia
using BayesMMfwd

configuration = BayesMMfwd.load_par_file("data/SIDES_from_original.par")
catalog = BayesMMfwd.load_sides_csv("data/SIDES_Bethermin2017_short2.csv")

inputs = BayesMMfwd.build_forward_inputs(catalog)
parameters = BayesMMfwd.build_forward_parameters(configuration; filters=false)
noise = BayesMMfwd.make_simulation_noise(length(inputs.redshift))

output = BayesMMfwd.forward_model(inputs, parameters.numeric, noise)
```

The explicit `noise` argument makes a run reproducible and lets all backends
consume identical random draws. The result is a named tuple of arrays; use
`BayesMMfwd.add_output_columns!` to attach those arrays to a catalog.

The Reactant path changes data placement and compilation, not the model:

```julia
using Reactant

reactant_inputs = Reactant.to_rarray(inputs)
reactant_parameters = Reactant.to_rarray(parameters.numeric)
reactant_noise = Reactant.to_rarray(noise)

compiled_forward = Reactant.@compile sync=true BayesMMfwd.forward_model(
    reactant_inputs,
    reactant_parameters,
    reactant_noise,
)
reactant_output = compiled_forward(
    reactant_inputs,
    reactant_parameters,
    reactant_noise,
)
```

## Production driver

`gen_pysides_from_original.jl` is a small Reactant driver around the forward
model. CPU is the default backend:

```bash
julia --project=. gen_pysides_from_original.jl
julia --project=. gen_pysides_from_original.jl --backend CUDA --large
julia --project=. gen_pysides_from_original.jl --backend TPU --large
```

It reports transfer, compilation, execution, materialization, and total time.
Add `--write-output` to write the configured FITS catalog, `--rows N` to
limit a run, or `--filters` to include configured filter fluxes. Run with
`--help` for all options.

## Performance and reproducibility

The Julia-versus-Reactant harness is isolated under [`benchmark/`](benchmark/).
It records machine-readable timing and allocation data and includes a Julia
plotting script for SVG and PDF figures.

## Tests

```bash
julia --startup-file=no --history-file=no --project=. test/runtests.jl
```

The regression suite checks Julia/Reactant agreement for the complete forward
model, input-table boundary cases, output assembly, and FITS writing.

## Source layout

- `src/reactant_pipeline.jl`: backend-independent numerical forward model.
- `src/pipeline_host.jl`: input preparation and table loading.
- `src/io.jl`: parameter, catalog, and FITS I/O.
- `data/`: tracked model inputs and local, ignored catalogs.
- `gen_pysides_from_original.jl`: Reactant execution and optional FITS output.
- `benchmark/`: benchmark runner and result utilities.
- `benchmark/plots/`: plotting scripts.
