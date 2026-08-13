using BayesMM: BayesMM
using Reactant: Reactant

include("common.jl")

function materialize_output(value::NamedTuple)
    names = keys(value)
    values = map(materialize_output, Tuple(value))
    return NamedTuple{names}(values)
end

materialize_output(value::Tuple) = map(materialize_output, value)
materialize_output(value::AbstractArray) = Array(value)
materialize_output(value) = value

function check_forward_output(reference, candidate; rtol=2.0e-6, atol=1.0e-8)
    keys(reference) == keys(candidate) ||
        error("Julia and Reactant returned different output fields")

    for name in keys(reference)
        expected = getproperty(reference, name)
        actual = getproperty(candidate, name)
        if expected === nothing
            actual === nothing || error("$name differs between Julia and Reactant")
            continue
        end

        size(expected) == size(actual) || error("$name has mismatched dimensions")
        if eltype(expected) <: Bool
            expected == actual || error("$name differs between Julia and Reactant")
            continue
        end

        all(isfinite, actual) || error("$name contains non-finite Reactant values")
        isapprox(expected, actual; rtol, atol) ||
            error("$name differs between Julia and Reactant")
    end

    @info "Correctness check passed" rtol atol
    return nothing
end

function run_bayesmm_benchmark!(
    results::Dict,
    backend::String;
    dataset="data/SIDES_Bethermin2017_short2.csv",
    param_path="data/SIDES_from_original.par",
    nrows=nothing,
    filters=false,
)
    params = BayesMM.load_par_file(param_path)
    catalog = BayesMM.load_sides_csv(dataset, nrows)
    inputs = BayesMM.build_forward_inputs(catalog)
    parameter_data = BayesMM.build_forward_parameters(params; filters)
    noise = BayesMM.make_simulation_noise(length(inputs.redshift))

    cpu_args = (inputs, parameter_data.numeric, noise)
    ra_args = (
        Reactant.to_rarray(inputs),
        Reactant.to_rarray(parameter_data.numeric),
        Reactant.to_rarray(noise),
    )

    n_gal = length(inputs.redshift)
    mode = filters ? "forward_filters" : "forward"
    benchmark_name = "BayesMM [$n_gal galaxies]/$mode"
    reactant_output = run_benchmark!(
        results,
        backend,
        benchmark_name,
        BayesMM.forward_model,
        cpu_args,
        ra_args;
        configs=[BenchmarkConfiguration("Default")],
    )
    reference = BayesMM.forward_model(cpu_args...)
    check_forward_output(reference, materialize_output(reactant_output))
    return n_gal
end
