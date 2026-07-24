include("gen_pysides_from_original.jl")

function run_julia_smoke(;
    dataset="data/SIDES_Bethermin2017_short2.csv",
    param_path="SIDES_from_original.par",
    nrows=nothing,
    filters=false,
    write_output=true,
)
    params = load_params(param_path)
    catalog = load_sides_csv(dataset, nrows)
    inputs = build_forward_inputs(catalog)
    parameter_data = build_forward_parameters(params; filters)
    noise = make_simulation_noise(length(inputs.redshift); seed=BENCHMARK_SEED)
    output = forward_model(inputs, parameter_data.numeric, noise)
    add_output_columns!(
        catalog,
        output,
        parameter_data.numeric.dust.lambda_list,
        parameter_data.filter_names,
    )
    write_output && gen_outputs(catalog, params)
    return catalog
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_julia_smoke()
end
