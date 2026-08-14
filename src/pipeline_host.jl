using CSV
using Cosmology
using DataFrames
using DelimitedFiles
using HDF5
using Random
using Unitful
using UnitfulAstro

const BENCHMARK_SEED = 1234
const SIDES_COSMOLOGY = Cosmology.cosmology(h=0.6774, OmegaM=0.3089)

float_param(params, key) = Float64(params[key])

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
