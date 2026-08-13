using InteractiveUtils: versioninfo

include("forward_model.jl")

function print_benchmark_usage(io=stdout)
    println(
        io,
        """
Usage:
  julia --project=benchmark benchmark/benchmark_reactant.jl [options]

Options:
  --dataset PATH       SIDES CSV input (default: data/SIDES_Bethermin2017_short2.csv)
  --rows N             Read only the first N rows
  --large              Use the large CSV (all rows unless followed by --rows)
  --filters            Include filter-grid fluxes (requires every FILTER_NAME.h5)
  --backend NAME       Reactant backend: CPU, CUDA, or TPU
                       (default: BENCHMARK_GROUP, then Reactant auto-detection)
  --results-dir PATH   JSON output directory (default: benchmark/results)
  --help               Show this message

Results use Reactant.jl's Runtime (s) and TFLOP/s schema. CPU runtime uses
Chairmarks for both Julia and Reactant; accelerator runtime and FLOP/s use
XProf. A separate memory JSON records host allocations and XProf backend-memory
peaks. Ordinary Julia is recorded only for CPU runs; every backend is checked
against it.
"""
    )
end

function parse_benchmark_command_line(args)
    options = (
        dataset="data/SIDES_Bethermin2017_short2.csv",
        nrows=nothing,
        filters=false,
        backend=nothing,
        results_dir=joinpath(@__DIR__, "results"),
    )
    values = Dict{Symbol,Any}(pairs(options))

    index = 1
    while index <= length(args)
        argument = args[index]
        if argument == "--help"
            print_benchmark_usage()
            return nothing
        elseif argument == "--large"
            values[:dataset] = "data/SIDES_Bethermin2017_short.csv"
            values[:nrows] = nothing
        elseif argument == "--filters"
            values[:filters] = true
        elseif argument in ("--dataset", "--rows", "--backend", "--results-dir")
            index == length(args) &&
                throw(ArgumentError("$argument requires a value"))
            value = args[index+1]
            if argument == "--dataset"
                values[:dataset] = value
            elseif argument == "--rows"
                values[:nrows] = parse(Int, value)
                values[:nrows] > 0 ||
                    throw(ArgumentError("--rows must be positive"))
            elseif argument == "--backend"
                group = uppercase(value)
                group == "GPU" && (group = "CUDA")
                group in ("CPU", "CUDA", "TPU") || throw(ArgumentError(
                    "--backend must be CPU, CUDA, or TPU",
                ))
                values[:backend] = group
            else
                values[:results_dir] = value
            end
            index += 1
        else
            throw(ArgumentError("unknown argument: $argument"))
        end
        index += 1
    end
    return (; values...)
end

function run_all_benchmarks(options)
    if !isnothing(options.backend)
        ENV["BENCHMARK_GROUP"] = options.backend
    end
    backend = get_backend()
    @info sprint(io -> versioninfo(io; verbose=true))

    results = Dict{String,Dict{String,Float64}}()
    n_gal = run_bayesmm_benchmark!(
        results,
        backend;
        dataset=options.dataset,
        nrows=options.nrows,
        filters=options.filters,
    )

    mode = options.filters ? "forward_filters" : "forward"
    prefix = "bayesmm_$(mode)_$n_gal"
    save_results(results, options.results_dir, prefix, backend)
    pretty_print_results(results, "BayesMM", backend)
    return results
end

function benchmark_main(args=ARGS)
    options = parse_benchmark_command_line(args)
    isnothing(options) && return nothing
    return run_all_benchmarks(options)
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_main()
end
