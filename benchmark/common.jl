using Reactant: Reactant
using Chairmarks: @b
using Printf: @sprintf

include("utils.jl")

struct BenchmarkConfiguration
    name::String
    compile_options::Union{Nothing,Reactant.CompileOptions}
    nrepeat::Int
end

function BenchmarkConfiguration(
    name::String;
    compile_options::Union{Nothing,Reactant.CompileOptions}=nothing,
    nrepeat::Int=25,
)
    return BenchmarkConfiguration(name, compile_options, nrepeat)
end

function run_benchmark!(
    results::Dict,
    backend::String,
    benchmark_name::String,
    fn::F,
    cpu_args::Tuple,
    ra_args::Tuple;
    skip_cpu::Bool=false,
    configs::Vector{BenchmarkConfiguration},
    benchmark_seconds::Float64=5.0,
    benchmark_samples::Int=100,
) where {F}
    for metric in (
        "Runtime (s)",
        "TFLOP/s",
        "Host allocated bytes",
        "Host allocations",
        "Backend peak memory (bytes)",
    )
        if !haskey(results, metric)
            results[metric] = Dict{String,Float64}()
        end
    end

    # Run CPU/Julia benchmark if on CPU backend
    if backend == "CPU" && !skip_cpu
        full_benchmark_name = string(benchmark_name, "/CPU/Julia")

        # Warmup
        fn(cpu_args...)

        # Benchmark using Chairmarks
        bench = @b fn($(cpu_args)...) seconds = benchmark_seconds evals = 1 samples =
            benchmark_samples

        results["Runtime (s)"][full_benchmark_name] = bench.time
        results["TFLOP/s"][full_benchmark_name] = -1.0
        results["Host allocated bytes"][full_benchmark_name] = Float64(bench.bytes)
        results["Host allocations"][full_benchmark_name] = Float64(bench.allocs)

        print_stmt = @sprintf(
            "%100s     :     %.5gs     %.5g MiB     %.0f allocations",
            full_benchmark_name,
            bench.time,
            bench.bytes / 2.0^20,
            bench.allocs,
        )
        @info print_stmt
        GC.gc(true)
    end

    # Run Reactant for each configuration. Chairmarks gives a directly
    # comparable CPU measurement; accelerators use synchronized XProf timing.
    output = nothing
    for config in configs
        full_benchmark_name = string(benchmark_name, "/", backend, "/", config.name)

        compile_options = something(config.compile_options, Reactant.CompileOptions())
        compile_options =
            Reactant.__compile_options_with_updated_sync(compile_options, true)

        if backend == "CPU"
            compiled_fn = Reactant.compile(fn, ra_args; compile_options)
            output = compiled_fn(ra_args...)
            allocation_bench = @b $compiled_fn($(ra_args)...) seconds =
                benchmark_seconds evals = 1 samples = benchmark_samples
            results["Runtime (s)"][full_benchmark_name] = allocation_bench.time
            results["TFLOP/s"][full_benchmark_name] = -1.0
        else
            prof_results = Reactant.Profiler.profile_with_xprof(
                fn,
                ra_args...;
                nrepeat=config.nrepeat,
                compile_options=config.compile_options,
            )
            output = prof_results.val
            results["Runtime (s)"][full_benchmark_name] =
                prof_results.profiling_result.runtime_ns / 1e9
            metrics = prof_results.profiling_result.metrics_data
            tflops = isnothing(metrics) ? -1.0 : metrics.raw_flops_rate / 1e12
            results["TFLOP/s"][full_benchmark_name] =
                isfinite(tflops) ? tflops : -1.0

            memory_data = prof_results.profiling_result.memory_data
            for allocator in sort!(collect(keys(memory_data)))
                memory_name = string(full_benchmark_name, " [", allocator, "]")
                peak_bytes = memory_data[allocator].peak_bytes_usage_lifetime
                results["Backend peak memory (bytes)"][memory_name] =
                    Float64(peak_bytes)
                @info "Reactant backend memory peak" benchmark =
                    full_benchmark_name allocator peak_bytes
            end

            output = nothing
            prof_results = nothing
            GC.gc(true)
            compiled_fn = Reactant.compile(fn, ra_args; compile_options)
            output = compiled_fn(ra_args...)
            allocation_bench = @b $compiled_fn($(ra_args)...) seconds =
                benchmark_seconds evals = 1 samples = 1
        end

        results["Host allocated bytes"][full_benchmark_name] =
            Float64(allocation_bench.bytes)
        results["Host allocations"][full_benchmark_name] =
            Float64(allocation_bench.allocs)

        print_stmt = @sprintf(
            "%100s     :     %.5gs     %.5g TFLOP/s     %.5g MiB     %.0f allocations",
            full_benchmark_name,
            results["Runtime (s)"][full_benchmark_name],
            results["TFLOP/s"][full_benchmark_name],
            allocation_bench.bytes / 2.0^20,
            allocation_bench.allocs,
        )
        @info print_stmt
        GC.gc(true)
    end

    return output
end
