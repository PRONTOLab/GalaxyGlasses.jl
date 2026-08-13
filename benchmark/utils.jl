# Adapted from Reactant.jl/benchmark/utils.jl so this benchmark can move into
# Reactant's benchmark suite without changing its result format.

using JSON: JSON
using PrettyTables: pretty_table

function standardized_results(values::Dict{String,Float64}, unit::String)
    results = Vector{Dict{String,Union{String,Float64}}}()
    for name in sort!(collect(keys(values)))
        push!(results, Dict("name" => name, "value" => values[name], "unit" => unit))
    end
    return results
end

function get_backend()
    benchmark_group = get(ENV, "BENCHMARK_GROUP", nothing)

    if benchmark_group == "CUDA"
        Reactant.set_default_backend("gpu")
        @info "Running CUDA benchmarks" maxlog = 1
    elseif benchmark_group == "TPU"
        Reactant.set_default_backend("tpu")
        @info "Running TPU benchmarks" maxlog = 1
    elseif benchmark_group == "CPU"
        Reactant.set_default_backend("cpu")
        @info "Running CPU benchmarks" maxlog = 1
    else
        benchmark_group = String(split(string(first(Reactant.devices())), ":")[1])
        @info "Running $(benchmark_group) benchmarks" maxlog = 1
    end

    @assert benchmark_group in ("CPU", "CUDA", "TPU") "Unknown backend: $(benchmark_group)"
    return benchmark_group
end

function save_results(
    results::Dict{String,Dict{String,Float64}},
    results_dir::String,
    prefix::String,
    backend::String,
)
    mkpath(results_dir)

    benchmark_filename = string(prefix, "_", backend, "benchmarks.json")
    tflops_filename = string(prefix, "_", backend, "benchmarks_tflops.json")
    memory_filename = string(prefix, "_", backend, "benchmarks_memory.json")
    benchmark_filepath = joinpath(results_dir, benchmark_filename)
    tflops_filepath = joinpath(results_dir, tflops_filename)
    memory_filepath = joinpath(results_dir, memory_filename)

    runtime_results = standardized_results(results["Runtime (s)"], "s")
    tflops_results = standardized_results(results["TFLOP/s"], "TFLOP/s")
    memory_results = Vector{Dict{String,Union{String,Float64}}}()
    memory_metrics = (
        ("Host allocated bytes", "B"),
        ("Host allocations", "allocations"),
        ("Backend peak memory (bytes)", "B"),
    )
    for (metric, unit) in memory_metrics
        values = get(results, metric, Dict{String,Float64}())
        for name in sort!(collect(keys(values)))
            push!(
                memory_results,
                Dict(
                    "name" => name,
                    "metric" => metric,
                    "value" => values[name],
                    "unit" => unit,
                ),
            )
        end
    end

    open(benchmark_filepath, "w") do io
        JSON.json(io, runtime_results; pretty=true)
    end
    open(tflops_filepath, "w") do io
        JSON.json(io, tflops_results; pretty=true)
    end
    open(memory_filepath, "w") do io
        JSON.json(io, memory_results; pretty=true)
    end

    @info "Saved benchmark results" runtime = benchmark_filepath tflops =
        tflops_filepath memory = memory_filepath
    return benchmark_filepath, tflops_filepath, memory_filepath
end

function pretty_print_results(results::Dict, suite::String, backend::String)
    runtime_results = get(results, "Runtime (s)", Dict{String,Float64}())
    tflops_results = get(results, "TFLOP/s", Dict{String,Float64}())
    byte_results = get(results, "Host allocated bytes", Dict{String,Float64}())
    allocation_results = get(results, "Host allocations", Dict{String,Float64}())

    if isempty(runtime_results)
        @warn "No benchmark results to display for $(suite)/$(backend)"
        return nothing
    end

    sorted_keys = sort(collect(keys(runtime_results)))
    table = Matrix{Any}(undef, length(sorted_keys), 8)
    for (index, key) in enumerate(sorted_keys)
        parts = split(key, "/")
        while length(parts) < 4
            push!(parts, "")
        end
        table[index, 1] = parts[1]
        table[index, 2] = parts[2]
        table[index, 3] = parts[3]
        table[index, 4] = join(parts[4:end], "/")
        table[index, 5] = runtime_results[key]
        table[index, 6] = get(tflops_results, key, -1.0)
        table[index, 7] = round(Int, get(byte_results, key, 0.0))
        table[index, 8] = round(Int, get(allocation_results, key, 0.0))
    end

    println()
    println("="^120)
    println("  Benchmark Results: $(suite) / $(backend)")
    println("="^120)
    pretty_table(
        table;
        alignment=[:l, :l, :l, :l, :c, :c, :c, :c],
        column_labels=[
            "Benchmark",
            "Mode",
            "Backend",
            "Config",
            "Time (s)",
            "TFLOP/s",
            "Host allocated (B)",
            "Host allocations",
        ],
        display_size=(-1, -1),
    )
    println("="^120)
    return nothing
end
