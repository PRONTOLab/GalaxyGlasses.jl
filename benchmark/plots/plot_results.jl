#!/usr/bin/env julia

using CairoMakie
using JSON
using Printf

# Keep pre-rename result files readable.
const RESULT_FILE =
    r"^bayesmm(?:fwd)?_(?:forward|forward_filters)_(\d+)_(CPU|CUDA|TPU)benchmarks\.json$"
const IMPLEMENTATIONS = ("Julia", "Reactant CPU", "Reactant CUDA", "Reactant TPU")
const LABELS = Dict(
    "Julia" => "Julia",
    "Reactant CPU" => "Reactant CPU",
    "Reactant CUDA" => "Reactant GPU",
    "Reactant TPU" => "Reactant TPU",
)
const COLORS = Dict(
    "Julia" => RGBf(0.298, 0.471, 0.659),
    "Reactant CPU" => RGBf(0.961, 0.522, 0.094),
    "Reactant CUDA" => RGBf(0.329, 0.635, 0.294),
    "Reactant TPU" => RGBf(0.698, 0.475, 0.635),
)
const MARKERS = Dict(
    "Julia" => :circle,
    "Reactant CPU" => :rect,
    "Reactant CUDA" => :utriangle,
    "Reactant TPU" => :diamond,
)
const LINESTYLES = Dict(
    "Julia" => :solid,
    "Reactant CPU" => :dash,
    "Reactant CUDA" => :dot,
    "Reactant TPU" => :dashdot,
)

function implementation(name)
    parts = split(name, '/')
    length(parts) >= 4 || error("unexpected benchmark name: $name")
    backend, configuration = parts[(end - 1):end]
    return configuration == "Julia" ? "Julia" : "Reactant $backend"
end

function add_measurement!(measurements, key, value)
    push!(get!(measurements, key, Float64[]), Float64(value))
end

function read_measurements(results_dir)
    runtimes = Dict{Tuple{Int,String},Vector{Float64}}()
    host_bytes = Dict{Tuple{Int,String},Vector{Float64}}()
    backend_peaks = Dict{Tuple{Int,String},Vector{Float64}}()

    for (root, _, files) in walkdir(results_dir)
        for filename in sort(files)
            matched = match(RESULT_FILE, filename)
            isnothing(matched) && continue

            rows = parse(Int, matched.captures[1])
            path = joinpath(root, filename)
            for result in JSON.parsefile(path)
                key = (rows, implementation(result["name"]))
                add_measurement!(runtimes, key, result["value"])
            end

            memory_path = replace(path, "benchmarks.json" => "benchmarks_memory.json")
            isfile(memory_path) || continue
            for result in JSON.parsefile(memory_path)
                name = replace(result["name"], r" \[[^][]+\]$" => "")
                key = (rows, implementation(name))
                metric = result["metric"]
                if metric == "Host allocated bytes"
                    add_measurement!(host_bytes, key, result["value"])
                elseif metric == "Backend peak memory (bytes)"
                    add_measurement!(backend_peaks, key, result["value"])
                end
            end
        end
    end

    isempty(runtimes) && error("no BayesMMfwd runtime JSON files found below $results_dir")
    return (; runtimes, host_bytes, backend_peaks)
end

function sample_median(values)
    sorted = sort(values)
    middle = length(sorted) ÷ 2
    return isodd(length(sorted)) ? sorted[middle + 1] : (sorted[middle] + sorted[middle + 1]) / 2
end

aggregate(measurements) = Dict(key => sample_median(values) for (key, values) in measurements)

function points(values, name)
    selected = sort(
        [(rows, value) for ((rows, implementation_name), value) in values if
         implementation_name == name && value > 0];
        by=first,
    )
    return first.(selected), last.(selected)
end

function records(measurements, field)
    medians = aggregate(measurements)
    record_type = NamedTuple{(:galaxies, :implementation, field, :process_runs)}
    return [
        record_type((rows, get(LABELS, name, name), value, length(measurements[(rows, name)]))) for
        ((rows, name), value) in sort(collect(medians); by=first)
    ]
end

function write_summary(path, results_dir, measurements)
    runtime_medians = aggregate(measurements.runtimes)
    julia = Dict(rows => value for ((rows, name), value) in runtime_medians if name == "Julia")
    speedups = [
        (
            galaxies=rows,
            implementation=get(LABELS, name, name),
            speedup_vs_julia=julia[rows] / value,
        ) for ((rows, name), value) in sort(collect(runtime_medians); by=first) if
        name != "Julia" && haskey(julia, rows) && value > 0
    ]
    summary = (
        results_directory=results_dir,
        runtime=records(measurements.runtimes, :median_seconds),
        speedup=speedups,
        host_allocations=records(measurements.host_bytes, :median_allocated_bytes),
        backend_peak_memory=records(measurements.backend_peaks, :median_peak_bytes),
    )

    mkpath(dirname(path))
    open(path, "w") do io
        JSON.print(io, summary, 2)
        println(io)
    end
end

function format_rows(value)
    value >= 1_000_000 && return @sprintf("%.3gM", value / 1_000_000)
    value >= 1_000 && return @sprintf("%.3gk", value / 1_000)
    return string(round(Int, value))
end

function format_seconds(value)
    value >= 1 && return @sprintf("%.3g s", value)
    value >= 1e-3 && return @sprintf("%.3g ms", value * 1e3)
    return @sprintf("%.3g μs", value * 1e6)
end

function format_bytes(value)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    unit = 1
    while value >= 1024 && unit < length(units)
        value /= 1024
        unit += 1
    end
    return @sprintf("%.3g %s", value, units[unit])
end

function add_series!(axis, values, names; legend=false)
    for name in names
        rows, samples = points(values, name)
        isempty(rows) && continue
        options = (
            color=COLORS[name],
            linestyle=LINESTYLES[name],
            marker=MARKERS[name],
            linewidth=2.5,
            markersize=10,
        )
        if legend
            scatterlines!(axis, rows, samples; options..., label=LABELS[name])
        else
            scatterlines!(axis, rows, samples; options...)
        end
    end
end

function make_plot(svg_path, pdf_path, measurements)
    runtimes = aggregate(measurements.runtimes)
    host_bytes = aggregate(measurements.host_bytes)
    names = [name for name in IMPLEMENTATIONS if any(key[2] == name for key in keys(runtimes))]
    all_rows = sort(unique(first.(keys(runtimes))))
    julia = Dict(rows => value for ((rows, name), value) in runtimes if name == "Julia")
    speedups = Dict(
        (rows, name) => julia[rows] / value for ((rows, name), value) in runtimes if
        name != "Julia" && haskey(julia, rows) && value > 0
    )

    figure = Figure(size=(1500, 500), fontsize=13)
    header = figure[1, 1:3] = GridLayout()
    Label(header[1, 1], "BayesMMfwd forward-model scaling"; fontsize=22, font=:bold)
    common = (
        xscale=log10,
        yscale=log10,
        xticks=(Float64.(all_rows), format_rows.(all_rows)),
        xgridcolor=(:gray, 0.25),
        ygridcolor=(:gray, 0.25),
    )
    runtime_axis = Axis(
        figure[2, 1];
        common...,
        title="Runtime (lower is better)",
        xlabel="Catalog size (galaxies)",
        ylabel="Model runtime",
        ytickformat=values -> format_seconds.(values),
    )
    speedup_axis = Axis(
        figure[2, 2];
        common...,
        title="Speedup over Julia (higher is better)",
        xlabel="Catalog size (galaxies)",
        ylabel="Speedup",
        ytickformat=values -> [@sprintf("%.3g×", value) for value in values],
    )
    allocation_axis = Axis(
        figure[2, 3];
        common...,
        title="Host allocations (lower is better)",
        xlabel="Catalog size (galaxies)",
        ylabel="Allocated bytes per call",
        ytickformat=values -> format_bytes.(values),
    )

    add_series!(runtime_axis, runtimes, names; legend=true)
    add_series!(speedup_axis, speedups, names)
    add_series!(allocation_axis, host_bytes, names)
    hlines!(speedup_axis, [1.0]; color=(:black, 0.45), linestyle=:dash, linewidth=1)
    Legend(header[2, 1], runtime_axis; orientation=:horizontal, framevisible=false)
    footer = string(
        "Warm execution only; compilation, transfers, and materialization excluded. ",
        "One process-level run per point.",
    )
    Label(
        figure[3, 1:3],
        footer;
        fontsize=11,
        color=:gray35,
    )
    colgap!(figure.layout, 24)
    rowgap!(figure.layout, 8)

    for path in (svg_path, pdf_path)
        mkpath(dirname(path))
        save(path, figure)
        println("wrote $path")
    end
end

function parse_options(args)
    options = Dict(
        "--results-dir" => "benchmark/results/paper-2026-08-13",
        "--svg" => "benchmark/bayesmmfwd_scaling.svg",
        "--pdf" => "benchmark/bayesmmfwd_scaling.pdf",
        "--summary" => "benchmark/results_summary.json",
    )
    if "--help" in args
        println("Usage: julia --project=benchmark benchmark/plots/plot_results.jl [options]")
        println("  --results-dir DIR   benchmark JSON directory")
        println("  --svg FILE          SVG output")
        println("  --pdf FILE          vector PDF output for LaTeX")
        println("  --summary FILE      aggregated JSON output")
        exit()
    end

    index = 1
    while index <= length(args)
        option = args[index]
        haskey(options, option) || error("unknown option: $option")
        index < length(args) || error("missing value for $option")
        options[option] = args[index + 1]
        index += 2
    end
    return options
end

options = parse_options(ARGS)
measurements = read_measurements(options["--results-dir"])
make_plot(options["--svg"], options["--pdf"], measurements)
write_summary(options["--summary"], options["--results-dir"], measurements)
println("wrote $(options["--summary"])")
