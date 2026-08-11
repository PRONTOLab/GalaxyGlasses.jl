using Base.Meta  # Needed for Meta.parse (optional, as Base.Meta exports Meta.parse)
using DataFrames
using CSV
using FITSIO
using Unitful
using UnitfulAstro
using Serialization

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

function gen_outputs(cat::DataFrame, params::Dict)
    # 1. Handle directory creation
    output_path = params["output_path"]
    if !isdir(output_path)
        println("Create $output_path")
        mkpath(output_path) # Equivalent to os.makedirs
    end

    # 2. Export to Pickle (Serialization)
    # In Julia, serialize is the standard way to dump objects to a binary file
    #if get(params, "gen_pickle", false) == true
    #    file_p = joinpath(output_path, params["run_name"] * ".p")
    #    println("Export the catalog to pickle... ($file_p)")
    #    open(file_p, "w") do f
    #        serialize(f, cat)
    #    end
    #end

    # 3. Export to FITS
    if get(params, "gen_fits", false) == true
        file_fits = joinpath(output_path, params["run_name"] * ".fits")
        println("Export the catalog to FITS... ($file_fits)")

        # Create a copy to add units without mutating the original catalog
        export_cat = copy(cat)
        col_names = names(export_cat)

        # Create a compatible dictionary for FITSIO
        data_dict = Dict{String,AbstractVector}()

        for name in col_names
            # 1. Extract the column and strip any physical units 
            col = ustrip.(export_cat[!, name])

            # 2. Check for BitVector and convert to standard Vector{Bool}
            # Standard Arrays are recognized by FITSIO's Array{T} method 
            if col isa BitVector
                col = Vector{Bool}(col)
            end

            data_dict[name] = col
        end

        # 1. Prepare three vectors to define the header records
        # FITS records consist of a Key, a Value, and a Comment
        h_keys = String[]
        h_vals = Any[]
        h_comms = String[]

        # 2. Populate the vectors with your simulation parameters
        for (key, val) in params
            # We use "COMMENT" as the keyword for every parameter
            push!(h_keys, "COMMENT")
            # COMMENT cards have no value in the standard FITS format
            push!(h_vals, nothing)
            # Store the "Key = Value" string in the comment field of the record
            push!(h_comms, "$key = $val")
        end
        header = FITSHeader(h_keys, h_vals, h_comms)

        # 2. Add your parameters as "COMMENT" cards to the header object
        # FITS supports multiple entries under the "COMMENT" keyword
        #for (key, val) in params
        #    # set_comment! adds or modifies metadata in a header object [1]
        #    # We use "COMMENT" as the key to create standard FITS comment lines
        #    set_comment!(header, "COMMENT", "$key = $val")
        #end

        # 3. Write to FITS using the function-block syntax to ensure the file closes
        FITS(file_fits, "w") do f
            # This will now succeed because all types are standard arrays
            write(f, data_dict; header=header)

            ## 4. Add simulation parameters as comments
            #for (key, val) in params
            #    FITSIO.write_key(f[1], "$key = $val")
            #end
        end
    end

    return true
end
