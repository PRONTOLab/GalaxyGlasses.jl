using Base.Meta  # Needed for Meta.parse (optional, as Base.Meta exports Meta.parse)
using DataFrames
using CSV

"""
    load_par_file(path::AbstractString, force_pysides_path::AbstractString="")

Reads a configuration file, ignoring comments starting with '#', and parses 
the values dynamically using Meta.parse and Base.eval, similar to Python's eval.
Returns a Dict{String, Any}.
"""
function load_par_file(path::AbstractString, force_pysides_path::AbstractString="")
    params = Dict{String,Any}()

    # 1. Read and Parse Key-Value Pairs
    # Use the 'do' block syntax for safe file handling (ensures 'file' is closed)
    open(path, "r") do file
        for line in eachline(file)
            line_stripped = strip(line)
            # Skip lines starting with '#' (comments)
            if !startswith(line_stripped, "#")
                # Split line by the first occurrence of '#', keeping only the key/value part
                no_comment = split(line_stripped, '#', limit=2)[1]
                # Split key and value by the first occurrence of '='
                key_value = split(no_comment, '=', limit=2)
                if length(key_value) == 2
                    key = strip(key_value[1])
                    value_str = strip(key_value[2])
                    params[key] = value_str
                end
            end
        end
    end

    # 2. Evaluate and Convert Types (Equivalent to Python's eval loop)
    # This dynamic evaluation is necessary to convert strings like "3.14" 
    # or "[1-3]" into their corresponding Julia types.
    for (key, value_str) in params
        try
            ex = Meta.parse(value_str)
            params[key] = Base.eval(Main, ex)
        catch e
            @warn "Could not evaluate parameter '$key'. Keeping value as raw string: $value_str" exception = (e, catch_backtrace())
        end
    end
    return params
end

"""
    load_sides_csv(catfile::AbstractString, nrows::Union{Nothing, Int} = nothing)

Loads a catalog CSV, selecting the primary columns required for analysis.
nrows: Limits the number of rows read from the file (excluding headers).
"""
function load_sides_csv(catfile::AbstractString, nrows::Union{Nothing,Int}=nothing)

    # Equivalent to Python's print()
    println("Load the catalog CSV generated from the original IDL code to get RA, Dec, z, Mhalo, and Mstar...") # [9]

    # Equivalent to pd.read_csv, reading directly into a DataFrame [8].
    # The 'limit' keyword argument handles the Python 'nrows' parameter [11].
    # We explicitly specify the delimiter as ',' using 'delim' (equivalent to 'sep=', assuming a standard CSV).
    cat_IDL = CSV.read(
        catfile,
        DataFrame;
        delim=',',
        limit=nrows # Reads only up to `nrows`
    )

    # --- Column Selection and Renaming (Equivalent to pd.DataFrame(..., columns=...)) ---

    # Define the list of required columns (as Symbols in Julia, which represent column names)
    required_cols = [:redshift, :ra, :dec, :Mhalo, :Mstar]

    # Select only the required columns and return a new DataFrame [2, 12].
    # This assumes the CSV file already contains columns with these exact names.
    cat = cat_IDL[!, required_cols]

    return cat
end
