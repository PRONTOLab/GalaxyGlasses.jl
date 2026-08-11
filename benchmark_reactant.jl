include("gen_pysides_from_original.jl")

function print_benchmark_usage(io=stdout)
    println(io, """
Usage:
  julia --project=. benchmark_reactant.jl [options]

Options:
  --dataset PATH       SIDES CSV input (default: data/SIDES_Bethermin2017_short2.csv)
  --rows N             Read only the first N rows
  --large              Use the large CSV (all rows unless followed by --rows)
  --filters            Include filter-grid fluxes (requires every FILTER_NAME.h5)
  --samples N          Number of synchronized warm samples (default: 5)
  --results PATH       Human-readable report (default: benchmark_results_reactant.txt)
  --csv-results PATH   Machine-readable row (default: same name with .csv)
  --help               Show this message

The benchmark gives Julia and Reactant the same inputs and random draws,
synchronizes Reactant before stopping each timer, and checks every output.
""")
end

function parse_benchmark_command_line(args)
    options = (
        dataset="data/SIDES_Bethermin2017_short2.csv",
        nrows=nothing,
        filters=false,
        samples=5,
        results_path="benchmark_results_reactant.txt",
        csv_results_path=nothing,
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
        elseif argument in (
            "--dataset",
            "--rows",
            "--samples",
            "--results",
            "--csv-results",
        )
            index == length(args) &&
                throw(ArgumentError("$argument requires a value"))
            value = args[index + 1]
            if argument == "--dataset"
                values[:dataset] = value
            elseif argument == "--rows"
                values[:nrows] = parse(Int, value)
                values[:nrows] > 0 ||
                    throw(ArgumentError("--rows must be positive"))
            elseif argument == "--samples"
                values[:samples] = parse(Int, value)
                values[:samples] > 0 ||
                    throw(ArgumentError("--samples must be positive"))
            elseif argument == "--results"
                values[:results_path] = value
            else
                values[:csv_results_path] = value
            end
            index += 1
        else
            throw(ArgumentError("unknown argument: $argument"))
        end
        index += 1
    end

    if isnothing(values[:csv_results_path])
        base, _ = splitext(values[:results_path])
        values[:csv_results_path] = base * ".csv"
    end
    return (; values...)
end

function run_benchmark(;
    dataset="data/SIDES_Bethermin2017_short2.csv",
    param_path="SIDES_from_original.par",
    nrows=nothing,
    filters=false,
    samples=5,
    results_path="benchmark_results_reactant.txt",
    csv_results_path="benchmark_results_reactant.csv",
)
    params = load_par_file(param_path)
    catalog = load_sides_csv(dataset, nrows)
    inputs = build_forward_inputs(catalog)
    parameter_data = build_forward_parameters(params; filters)
    n_gal = length(inputs.redshift)
    noise = make_simulation_noise(n_gal; seed=BENCHMARK_SEED)

    println("Benchmarking Julia and Reactant for $n_gal rows...")
    benchmark = benchmark_forward_model(
        inputs,
        parameter_data.numeric,
        noise;
        samples,
    )
    filter_count = length(parameter_data.filter_names)
    report = format_benchmark_report(
        dataset,
        n_gal,
        benchmark.timings,
        benchmark.correctness;
        samples,
        filter_count,
    )
    print(report)
    write(results_path, report)

    row = benchmark_result_row(
        dataset,
        n_gal,
        benchmark.timings,
        benchmark.correctness;
        samples,
        filter_count,
    )
    CSV.write(csv_results_path, DataFrame([row]))
    println("Text report: $results_path")
    println("CSV result:  $csv_results_path")
    return (; benchmark, report, row)
end

function benchmark_main(args=ARGS)
    options = parse_benchmark_command_line(args)
    isnothing(options) && return nothing
    return run_benchmark(; options...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_main()
end
