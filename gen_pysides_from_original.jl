using Reactant
using DataFrames

Reactant.set_default_backend("cpu")

include("load_params.jl")
include("load_sides_csv.jl")
include("reactant_pipeline.jl")
include("pipeline_host.jl")
include("gen_outputs.jl")

function print_usage(io=stdout)
    println(io, """
Usage:
  julia --project=. gen_pysides_from_original.jl [options]

Options:
  --dataset PATH     SIDES CSV input (default: data/SIDES_Bethermin2017_short2.csv)
  --rows N           Read only the first N rows
  --large            Use the large CSV (all rows unless followed by --rows)
  --filters          Include filter-grid fluxes (requires every FILTER_NAME.h5)
  --write-output     Write the Reactant result using gen_outputs.jl
  --help             Show this message

This is the Reactant-only production driver. For Julia-versus-Reactant
comparisons, use benchmark_reactant.jl.
""")
end

function parse_command_line(args)
    options = (
        dataset="data/SIDES_Bethermin2017_short2.csv",
        nrows=nothing,
        filters=false,
        write_output=false,
    )
    values = Dict{Symbol,Any}(pairs(options))

    index = 1
    while index <= length(args)
        argument = args[index]
        if argument == "--help"
            print_usage()
            return nothing
        elseif argument == "--large"
            values[:dataset] = "data/SIDES_Bethermin2017_short.csv"
            values[:nrows] = nothing
        elseif argument == "--filters"
            values[:filters] = true
        elseif argument == "--write-output"
            values[:write_output] = true
        elseif argument in ("--dataset", "--rows")
            index == length(args) &&
                throw(ArgumentError("$argument requires a value"))
            value = args[index + 1]
            if argument == "--dataset"
                values[:dataset] = value
            else
                values[:nrows] = parse(Int, value)
                values[:nrows] > 0 ||
                    throw(ArgumentError("--rows must be positive"))
            end
            index += 1
        else
            throw(ArgumentError("unknown argument: $argument"))
        end
        index += 1
    end
    return (; values...)
end

function run_reactant_pipeline(;
    dataset="data/SIDES_Bethermin2017_short2.csv",
    param_path="SIDES_from_original.par",
    nrows=nothing,
    filters=false,
    write_output=false,
)
    params = load_params(param_path)
    catalog_template = load_sides_csv(dataset, nrows)
    inputs = build_forward_inputs(catalog_template)
    parameter_data = build_forward_parameters(params; filters)
    noise = make_simulation_noise(length(inputs.redshift); seed=BENCHMARK_SEED)

    n_gal = length(inputs.redshift)
    println("Running the Reactant forward model for $n_gal rows...")
    reactant_run = run_reactant_forward_model(
        inputs,
        parameter_data.numeric,
        noise,
    )
    timings = reactant_run.timings
    println("Reactant timings:")
    println("  host-to-device: ", @sprintf("%.6f s", timings.host_to_device))
    println("  compilation:    ", @sprintf("%.6f s", timings.compilation))
    println("  execution:      ", @sprintf("%.6f s", timings.execution))
    println("  device-to-host: ", @sprintf("%.6f s", timings.device_to_host))
    println("  total:          ", @sprintf("%.6f s", timings.total))

    output_catalog = add_output_columns!(
        copy(catalog_template),
        reactant_run.output,
        parameter_data.numeric.dust.lambda_list,
        parameter_data.filter_names,
    )
    if write_output
        gen_outputs(output_catalog, params)
    end

    return (; output_catalog, reactant_run, params)
end

function main(args=ARGS)
    options = parse_command_line(args)
    isnothing(options) && return nothing
    return run_reactant_pipeline(; options...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
