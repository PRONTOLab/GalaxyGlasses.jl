import BayesMM

function run_julia_smoke(;
    dataset="data/SIDES_Bethermin2017_short2.csv",
    param_path="SIDES_from_original.par",
    nrows=nothing,
    filters=false,
    write_output=true,
)
    params = BayesMM.load_par_file(param_path)
    catalog = BayesMM.load_sides_csv(dataset, nrows)
    inputs = BayesMM.build_forward_inputs(catalog)
    parameter_data = BayesMM.build_forward_parameters(params; filters)
    noise = BayesMM.make_simulation_noise(length(inputs.redshift))
    output = BayesMM.forward_model(inputs, parameter_data.numeric, noise)
    BayesMM.add_output_columns!(
        catalog,
        output,
        parameter_data.numeric.dust.lambda_list,
        parameter_data.filter_names,
    )
    write_output && BayesMM.gen_outputs(catalog, params)
    return catalog
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_julia_smoke()
end
