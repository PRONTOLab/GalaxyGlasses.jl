using CSV
using Cosmology
using DataFrames
using Dates
using DelimitedFiles
using HDF5
using Printf
using Random
using Statistics
using Unitful
using UnitfulAstro

const BENCHMARK_SEED = 1234
const SIDES_COSMOLOGY = Cosmology.cosmology(h=0.6774, OmegaM=0.3089)

float_param(params, key) = Float64(params[key])
cpu_model() = isempty(Sys.cpu_info()) ? Sys.CPU_NAME : Sys.cpu_info()[1].model

function parse_sfr_numeric_params(params)
    return (
        Chab2Salp_num=float_param(params, "Chab2Salp"),
        Mt0=float_param(params, "Mt0"),
        alpha1=float_param(params, "alpha1"),
        alpha2=float_param(params, "alpha2"),
        sigma0=float_param(params, "sigma0"),
        beta1=float_param(params, "beta1"),
        beta2=float_param(params, "beta2"),
        qfrac0=float_param(params, "qfrac0"),
        gamma=float_param(params, "gamma"),
        m1=float_param(params, "m1"),
        a2=float_param(params, "a2"),
        m0=float_param(params, "m0"),
        a0=float_param(params, "a0"),
        a1=float_param(params, "a1"),
        corr_zmean_lowzcorr=float_param(params, "corr_zmean_lowzcorr"),
        zmax_lowzcorr=float_param(params, "zmax_lowzcorr"),
        zmean_lowzcorr=float_param(params, "zmean_lowzcorr"),
        Psb_hz=float_param(params, "Psb_hz"),
        slope_Psb=float_param(params, "slope_Psb"),
        z_Psb_knee=float_param(params, "z_Psb_knee"),
        sigma_MS=float_param(params, "sigma_MS"),
        logx0=float_param(params, "logx0"),
        logBsb=float_param(params, "logBsb"),
        SFR_max=float_param(params, "SFR_max"),
    )
end

function read_hdf5_dict(path)
    isfile(path) || throw(ArgumentError("required HDF5 file does not exist: $path"))
    return h5open(path, "r") do file
        Dict(String(key) => read(file[key]) for key in keys(file))
    end
end

function orient_sed_matrix(matrix, n_lambda, n_u, name)
    if size(matrix) == (n_lambda, n_u)
        return matrix
    elseif size(matrix) == (n_u, n_lambda)
        return permutedims(matrix)
    end
    throw(DimensionMismatch(
        "$name has size $(size(matrix)); expected ($n_lambda, $n_u) or ($n_u, $n_lambda)",
    ))
end

function load_dust_tables(params)
    sed_dict = read_hdf5_dict(params["SED_file"])
    sed_lambda = Float64.(vec(sed_dict["lambda"]))
    sed_umean = Float64.(vec(sed_dict["Umean"]))
    n_lambda = length(sed_lambda)
    n_u = length(sed_umean)
    sed_tables = (
        lambda=sed_lambda,
        dU=Float64(sed_dict["dU"]),
        Umean=sed_umean,
        nuLnu_MS=Float64.(orient_sed_matrix(
            sed_dict["nuLnu_MS_arr"],
            n_lambda,
            n_u,
            "nuLnu_MS_arr",
        )),
        nuLnu_SB=Float64.(orient_sed_matrix(
            sed_dict["nuLnu_SB_arr"],
            n_lambda,
            n_u,
            "nuLnu_SB_arr",
        )),
    )

    ratio_dict = read_hdf5_dict(params["ratios_file"])
    ratio_umean = Float64.(vec(ratio_dict["Umean"]))
    n_ratio = length(ratio_umean)
    ratio_ms = Float64.(vec(ratio_dict["LFIR_LIR_ratio_MS"]))
    ratio_sb = Float64.(vec(ratio_dict["LFIR_LIR_ratio_SB"]))
    length(ratio_ms) == n_ratio ||
        throw(DimensionMismatch("LFIR MS ratio grid does not match its Umean grid"))
    length(ratio_sb) == n_ratio ||
        throw(DimensionMismatch("LFIR SB ratio grid does not match its Umean grid"))
    ratio_tables = (
        Umean=ratio_umean,
        dU=Float64(ratio_dict["dU"]),
        LFIR_LIR_ratio_MS=ratio_ms,
        LFIR_LIR_ratio_SB=ratio_sb,
    )

    return (
        UmeanSB=float_param(params, "UmeanSB"),
        UmeanMSz0=float_param(params, "UmeanMSz0"),
        alphaMS=float_param(params, "alphaMS"),
        zlimMS=float_param(params, "zlimMS"),
        sigma_logUmean=float_param(params, "sigma_logUmean"),
        SFR2LIR=float_param(params, "SFR2LIR"),
        lambda_list=Float64.(params["lambda_list"]),
        sed_tables,
        ratio_tables,
    )
end

function load_magnification_tables(path)
    isfile(path) || throw(ArgumentError("magnification grid does not exist: $path"))
    data_frame = CSV.read(path, DataFrame; comment="#")
    data = Matrix{Float64}(data_frame)
    z_grid = vec(data[1, 2:end])
    mu_grid = reverse(vec(data[2:end, 1]))
    psupmu = reverse(data[2:end, 2:end]; dims=1)
    size(psupmu) == (length(mu_grid), length(z_grid)) ||
        throw(DimensionMismatch("magnification grid dimensions are inconsistent"))
    return (; z_grid, mu_grid, Psupmu=psupmu)
end

function orient_filter_matrix(matrix, n_u, n_z, name)
    if size(matrix) == (n_u, n_z)
        return matrix
    elseif size(matrix) == (n_z, n_u)
        return permutedims(matrix)
    end
    throw(DimensionMismatch(
        "$name has size $(size(matrix)); expected ($n_u, $n_z) or ($n_z, $n_u)",
    ))
end

function load_filter_table(path)
    dictionary = read_hdf5_dict(path)
    umean = Float64.(vec(dictionary["Umean"]))
    n_u = length(umean)
    n_z = Int(dictionary["Nlog1plusz"])
    return (
        dU=Float64(dictionary["dU"]),
        dlog1plusz=Float64(dictionary["dlog1plusz"]),
        Nlog1plusz=n_z,
        Umean=umean,
        Snu_LIR_MS=Float64.(orient_filter_matrix(
            dictionary["Snu_LIR_MS"],
            n_u,
            n_z,
            "Snu_LIR_MS",
        )),
        Snu_LIR_SB=Float64.(orient_filter_matrix(
            dictionary["Snu_LIR_SB"],
            n_u,
            n_z,
            "Snu_LIR_SB",
        )),
    )
end

function load_filter_tables(params; enabled=false)
    names = String.(params["filter_list"])
    enabled || return (names=String[], tables=())

    grid_path = String(params["grid_filter_path"])
    paths = [joinpath(grid_path, name * ".h5") for name in names]
    missing_paths = filter(!isfile, paths)
    isempty(missing_paths) || throw(ArgumentError(
        "filter grids were requested but are missing:\n" * join(missing_paths, "\n"),
    ))
    tables = Tuple(load_filter_table(path) for path in paths)
    return (; names, tables)
end

function parse_line_numeric_params(params)
    sled_path = String(params["SLED_filename"])
    isfile(sled_path) || throw(ArgumentError("SLED file does not exist: $sled_path"))
    sled = Matrix{Float64}(readdlm(sled_path; comments=true, comment_char='#'))
    size(sled, 2) >= 3 ||
        throw(DimensionMismatch("SLED file must have at least three columns"))
    idiffuse = vec(sled[:, 2])
    iclump = vec(sled[:, 3])
    length(idiffuse) >= 8 ||
        throw(DimensionMismatch("SLED file must contain transitions J=1 through J=8"))
    idiffuse ./= idiffuse[1]
    iclump ./= iclump[1]

    return (
        sigma_dex_CO10=float_param(params, "sigma_dex_CO10"),
        nu_CO=float_param(params, "nu_CO"),
        Idiffuse=idiffuse,
        Iclump=iclump,
        SLED_SB_Birkin=Bool(params["SLED_SB_Birkin"]),
        rJup1_Birkin=Float64[
            0.9,
            0.6,
            0.32,
            0.35,
            0.3,
            0.22,
            0.22 * (7 / 8) ^ 2,
        ],
        nu_CII=float_param(params, "nu_CII"),
        sigma_dex_CII=float_param(params, "sigma_dex_CII"),
        a_CI10=float_param(params, "a_CI10"),
        b_CI10=float_param(params, "b_CI10"),
        sigma_CI10=float_param(params, "sigma_CI10"),
        a_CI21=float_param(params, "a_CI21"),
        b_CI21=float_param(params, "b_CI21"),
        sigma_CI21=float_param(params, "sigma_CI21"),
        nu_CI10=float_param(params, "nu_CI10"),
        nu_CI21=float_param(params, "nu_CI21"),
    )
end

function make_simulation_noise(n_gal; seed=BENCHMARK_SEED)
    rng = MersenneTwister(seed)
    return (
        q_uniform=rand(rng, n_gal),
        sb_uniform=rand(rng, n_gal),
        sfr_uniform=rand(rng, n_gal),
        mu_uniform=rand(rng, n_gal),
        umean_normal=randn(rng, n_gal),
        co_normal=randn(rng, n_gal),
        cii_lagache_normal=randn(rng, n_gal),
        cii_delooze_normal=randn(rng, n_gal),
        ci10_normal=randn(rng, n_gal),
        ci21_normal=randn(rng, n_gal),
    )
end

function build_forward_inputs(cat)
    redshift = Float64.(cat.redshift)
    dlum_mpc = [
        ustrip(u"Mpc", Cosmology.luminosity_dist(SIDES_COSMOLOGY, z))
        for z in redshift
    ]
    return (
        redshift,
        Mstar=Float64.(cat.Mstar),
        dlum_mpc,
        dlum_m=dlum_mpc .* MPC_TO_M,
    )
end

function build_forward_parameters(params; filters=false)
    filter_data = load_filter_tables(params; enabled=filters)
    numeric = (
        sfr=parse_sfr_numeric_params(params),
        magnification=load_magnification_tables(params["path_mu_file"]),
        dust=load_dust_tables(params),
        filters=filter_data.tables,
        lines=parse_line_numeric_params(params),
    )
    return (; numeric, filter_names=filter_data.names)
end

function materialize_output(value::NamedTuple)
    names = keys(value)
    values = map(materialize_output, Tuple(value))
    return NamedTuple{names}(values)
end

materialize_output(value::Tuple) = map(materialize_output, value)
materialize_output(value::AbstractArray) = Array(value)
materialize_output(value) = value

function output_arrays!(result, prefix, value::NamedTuple)
    for name in keys(value)
        child_prefix = isempty(prefix) ? String(name) : "$prefix.$name"
        output_arrays!(result, child_prefix, getproperty(value, name))
    end
    return result
end

function output_arrays!(result, prefix, value::AbstractArray)
    push!(result, prefix => value)
    return result
end

output_arrays!(result, prefix, value) = result

function compare_forward_outputs(reference, candidate; rtol=2.0e-6, atol=1.0e-8)
    reference_arrays = output_arrays!(Pair{String,Any}[], "", reference)
    candidate_arrays = Dict(output_arrays!(Pair{String,Any}[], "", candidate))
    details = NamedTuple[]
    all_passed = true

    for (name, expected) in reference_arrays
        actual = candidate_arrays[name]
        same_shape = size(expected) == size(actual)
        finite = if eltype(actual) <: AbstractFloat
            all(isfinite, actual)
        else
            true
        end
        passed = if !same_shape
            false
        elseif eltype(expected) <: Bool
            expected == actual
        else
            isapprox(expected, actual; rtol, atol)
        end

        max_absolute = if isempty(expected) || !same_shape || eltype(expected) <: Bool
            passed ? 0.0 : Inf
        else
            maximum(abs.(Float64.(actual) .- Float64.(expected)))
        end
        scale = if isempty(expected) || eltype(expected) <: Bool
            1.0
        else
            max(maximum(abs, Float64.(expected)), eps(Float64))
        end
        max_relative = max_absolute / scale
        push!(details, (; name, passed=passed && finite, max_absolute, max_relative))
        all_passed &= passed && finite
    end
    return (; passed=all_passed, details)
end

function add_output_columns!(cat, output, lambda_list, filter_names)
    cat[!, :qflag] = output.qflag
    cat[!, :SFR] = output.SFR
    cat[!, :issb] = output.issb
    cat[!, :mu] = output.mu
    cat[!, :Dlum] = [
        ustrip(u"Mpc", Cosmology.luminosity_dist(SIDES_COSMOLOGY, z))
        for z in cat.redshift
    ]
    cat[!, :Umean] = output.Umean
    cat[!, :LIR] = output.LIR
    for (index, wavelength) in enumerate(lambda_list)
        cat[!, Symbol("S$(round(Int, wavelength))")] =
            output.monochromatic_fluxes[:, index]
    end
    cat[!, :LFIR] = output.LFIR
    for (index, name) in enumerate(filter_names)
        cat[!, Symbol("S$name")] = output.filter_fluxes[:, index]
    end

    cat[!, :LprimCO10] = output.LprimCO10
    cat[!, :ICO10] = output.ICO10
    for k in 1:7
        cat[!, Symbol("ICO$(k + 1)$k")] = output.co_fluxes[:, k]
    end
    cat[!, :LCII_Lagache] = output.LCII_Lagache
    cat[!, :ICII_Lagache] = output.ICII_Lagache
    cat[!, :LCII_de_Looze] = output.LCII_de_Looze
    cat[!, :ICII_de_Looze] = output.ICII_de_Looze
    cat[!, :ICI10] = output.ICI10
    cat[!, :ICI21] = output.ICI21
    return cat
end

function seconds_since(start_ns)
    return (time_ns() - start_ns) / 1.0e9
end

function timed_samples(f, samples)
    times = Vector{Float64}(undef, samples)
    result = nothing
    for index in eachindex(times)
        start_ns = time_ns()
        result = f()
        times[index] = seconds_since(start_ns)
    end
    return (; result, times, median=median(times), minimum=minimum(times))
end

"""
    run_reactant_forward_model(inputs, params, noise)

Compile and execute the numerical model with Reactant only. Host-to-device
transfer, compilation, synchronized execution, and device-to-host
materialization are timed separately so production runs do not need to execute
the ordinary Julia baseline.
"""
function run_reactant_forward_model(inputs, params, noise)
    total_start = time_ns()

    transfer_start = time_ns()
    inputs_reactant = Reactant.to_rarray(inputs)
    params_reactant = Reactant.to_rarray(params)
    noise_reactant = Reactant.to_rarray(noise)
    host_to_device = seconds_since(transfer_start)

    compile_start = time_ns()
    compiled_forward = @compile sync=true forward_model(
        inputs_reactant,
        params_reactant,
        noise_reactant,
    )
    compilation = seconds_since(compile_start)

    execution_start = time_ns()
    reactant_output = compiled_forward(
        inputs_reactant,
        params_reactant,
        noise_reactant,
    )
    Reactant.synchronize(reactant_output)
    execution = seconds_since(execution_start)

    materialization_start = time_ns()
    host_output = materialize_output(reactant_output)
    device_to_host = seconds_since(materialization_start)
    total = seconds_since(total_start)

    timings = (; host_to_device, compilation, execution, device_to_host, total)
    return (; output=host_output, timings)
end

function benchmark_forward_model(inputs, params, noise; samples=5)
    samples > 0 || throw(ArgumentError("samples must be positive"))

    julia_start = time_ns()
    julia_first_output = forward_model(inputs, params, noise)
    julia_first = seconds_since(julia_start)
    julia_warm = timed_samples(
        () -> forward_model(inputs, params, noise),
        samples,
    )

    transfer_start = time_ns()
    inputs_reactant = Reactant.to_rarray(inputs)
    params_reactant = Reactant.to_rarray(params)
    noise_reactant = Reactant.to_rarray(noise)
    reactant_transfer = seconds_since(transfer_start)

    compile_start = time_ns()
    compiled_forward = @compile sync=true forward_model(
        inputs_reactant,
        params_reactant,
        noise_reactant,
    )
    reactant_compile = seconds_since(compile_start)

    first_start = time_ns()
    reactant_first_output = compiled_forward(
        inputs_reactant,
        params_reactant,
        noise_reactant,
    )
    Reactant.synchronize(reactant_first_output)
    reactant_first = seconds_since(first_start)

    first_materialization_start = time_ns()
    reactant_host_output = materialize_output(reactant_first_output)
    reactant_first_device_to_host = seconds_since(first_materialization_start)

    reactant_warm_times = Vector{Float64}(undef, samples)
    reactant_device_to_host_times = Vector{Float64}(undef, samples)
    reactant_materialized_times = Vector{Float64}(undef, samples)
    for index in 1:samples
        sample_start = time_ns()
        reactant_warm_output = compiled_forward(
            inputs_reactant,
            params_reactant,
            noise_reactant,
        )
        Reactant.synchronize(reactant_warm_output)
        reactant_warm_times[index] = seconds_since(sample_start)

        materialization_start = time_ns()
        materialize_output(reactant_warm_output)
        reactant_device_to_host_times[index] = seconds_since(materialization_start)
        reactant_materialized_times[index] = seconds_since(sample_start)
    end

    correctness = compare_forward_outputs(julia_first_output, reactant_host_output)
    if !correctness.passed
        failures = filter(detail -> !detail.passed, correctness.details)
        summary = join(
            (
                "$(detail.name) (rel=$(detail.max_relative), abs=$(detail.max_absolute))"
                for detail in failures
            ),
            ", ",
        )
        error("Reactant output does not match the Julia output: $summary")
    end

    timings = (
        julia_first,
        julia_warm_median=julia_warm.median,
        julia_warm_minimum=julia_warm.minimum,
        reactant_transfer,
        reactant_compile,
        reactant_first,
        reactant_first_device_to_host,
        reactant_warm_median=median(reactant_warm_times),
        reactant_warm_minimum=minimum(reactant_warm_times),
        reactant_device_to_host_median=median(reactant_device_to_host_times),
        reactant_device_to_host_minimum=minimum(reactant_device_to_host_times),
        reactant_materialized_median=median(reactant_materialized_times),
        reactant_materialized_minimum=minimum(reactant_materialized_times),
    )
    return (
        julia_output=julia_first_output,
        reactant_output=reactant_host_output,
        timings,
        correctness,
    )
end

function format_benchmark_report(
    dataset_path,
    n_gal,
    timings,
    correctness;
    samples,
    filter_count,
)
    resident_speedup = timings.julia_warm_median / timings.reactant_warm_median
    materialized_speedup =
        timings.julia_warm_median / timings.reactant_materialized_median
    reactant_end_to_end = timings.reactant_compile + timings.reactant_first
    reactant_one_shot = (
        timings.reactant_transfer +
        timings.reactant_compile +
        timings.reactant_first +
        timings.reactant_first_device_to_host
    )
    worst = sort(correctness.details; by=detail -> detail.max_relative, rev=true)

    io = IOBuffer()
    println(io, "BayesMMfwd Julia vs Reactant Benchmark")
    println(io, "Generated: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println(io, "Dataset: ", dataset_path)
    println(io, "Rows: ", n_gal)
    println(io, "Julia: ", VERSION)
    println(io, "Reactant: ", Base.pkgversion(Reactant))
    println(io, "Backend: cpu")
    println(io, "CPU: ", cpu_model())
    println(io, "CPU target: ", Sys.CPU_NAME)
    println(io, "Logical CPU threads: ", Sys.CPU_THREADS)
    println(io, "Julia threads: ", Threads.nthreads())
    println(io, "System: ", Sys.KERNEL, " ", Sys.MACHINE)
    println(io, "Seed: ", BENCHMARK_SEED)
    println(io, "Warm samples: ", samples)
    println(io, "Filter grids: ", filter_count)
    println(io, "Line stages: CO, CII, CI")
    println(io)
    println(io, rpad("Metric", 42), "Time (s)")
    println(io, repeat("-", 58))
    println(io, rpad("Julia first call (includes JIT)", 42),
            @sprintf("%.6f", timings.julia_first))
    println(io, rpad("Julia warm median", 42),
            @sprintf("%.6f", timings.julia_warm_median))
    println(io, rpad("Julia warm minimum", 42),
            @sprintf("%.6f", timings.julia_warm_minimum))
    println(io, rpad("Reactant host-to-device transfer", 42),
            @sprintf("%.6f", timings.reactant_transfer))
    println(io, rpad("Reactant compilation", 42),
            @sprintf("%.6f", timings.reactant_compile))
    println(io, rpad("Reactant first synchronized execution", 42),
            @sprintf("%.6f", timings.reactant_first))
    println(io, rpad("Reactant first device-to-host", 42),
            @sprintf("%.6f", timings.reactant_first_device_to_host))
    println(io, rpad("Reactant warm synchronized median", 42),
            @sprintf("%.6f", timings.reactant_warm_median))
    println(io, rpad("Reactant warm synchronized minimum", 42),
            @sprintf("%.6f", timings.reactant_warm_minimum))
    println(io, rpad("Reactant warm device-to-host median", 42),
            @sprintf("%.6f", timings.reactant_device_to_host_median))
    println(io, rpad("Reactant warm materialized median", 42),
            @sprintf("%.6f", timings.reactant_materialized_median))
    println(io, rpad("Reactant warm materialized minimum", 42),
            @sprintf("%.6f", timings.reactant_materialized_minimum))
    println(io, rpad("Reactant compile + first execution", 42),
            @sprintf("%.6f", reactant_end_to_end))
    println(io, rpad("Reactant one-shot materialized total", 42),
            @sprintf("%.6f", reactant_one_shot))
    println(io)
    println(io, "Resident Julia / Reactant speedup: ",
            @sprintf("%.3fx", resident_speedup))
    println(io, "Materialized Julia / Reactant speedup: ",
            @sprintf("%.3fx", materialized_speedup))
    println(io, "Correctness: PASS (rtol=2e-6, atol=1e-8)")
    println(io, "Worst relative differences:")
    for detail in Iterators.take(worst, min(5, length(worst)))
        println(
            io,
            "  ",
            rpad(detail.name, 34),
            @sprintf("rel=%.3e abs=%.3e", detail.max_relative, detail.max_absolute),
        )
    end
    return String(take!(io))
end

function benchmark_result_row(
    dataset_path,
    n_gal,
    timings,
    correctness;
    samples,
    filter_count,
)
    max_relative = maximum(detail.max_relative for detail in correctness.details)
    return (
        generated=Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
        dataset=dataset_path,
        rows=n_gal,
        julia_version=string(VERSION),
        reactant_version=string(Base.pkgversion(Reactant)),
        backend="cpu",
        cpu=cpu_model(),
        cpu_target=Sys.CPU_NAME,
        logical_cpu_threads=Sys.CPU_THREADS,
        julia_threads=Threads.nthreads(),
        system="$(Sys.KERNEL) $(Sys.MACHINE)",
        seed=BENCHMARK_SEED,
        warm_samples=samples,
        filter_grids=filter_count,
        julia_first_s=timings.julia_first,
        julia_warm_median_s=timings.julia_warm_median,
        julia_warm_minimum_s=timings.julia_warm_minimum,
        reactant_host_to_device_s=timings.reactant_transfer,
        reactant_compile_s=timings.reactant_compile,
        reactant_first_s=timings.reactant_first,
        reactant_first_device_to_host_s=timings.reactant_first_device_to_host,
        reactant_warm_median_s=timings.reactant_warm_median,
        reactant_warm_minimum_s=timings.reactant_warm_minimum,
        reactant_device_to_host_median_s=timings.reactant_device_to_host_median,
        reactant_device_to_host_minimum_s=timings.reactant_device_to_host_minimum,
        reactant_materialized_median_s=timings.reactant_materialized_median,
        reactant_materialized_minimum_s=timings.reactant_materialized_minimum,
        reactant_compile_plus_first_s=(
            timings.reactant_compile + timings.reactant_first
        ),
        reactant_one_shot_materialized_s=(
            timings.reactant_transfer +
            timings.reactant_compile +
            timings.reactant_first +
            timings.reactant_first_device_to_host
        ),
        julia_over_reactant_resident_speedup=(
            timings.julia_warm_median / timings.reactant_warm_median
        ),
        julia_over_reactant_materialized_speedup=(
            timings.julia_warm_median / timings.reactant_materialized_median
        ),
        correctness_passed=correctness.passed,
        maximum_relative_difference=max_relative,
    )
end
