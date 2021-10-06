"""
This function returns a dataframe of given csv data.
"""
function read_data(file_name::String)
    projectdata = DataFrames.DataFrame(CSV.File(file_name;
                        truestrings=["T", "TRUE", "true"],
                        falsestrings=["F", "FALSE", "false"]));
    return projectdata
end

function write_data(dir:: String, file_name::String, data::DataFrames.DataFrame)

    dir_exists(dir::String)
    CSV.write(joinpath(dir, file_name), data)

    return
end

"""
This function makes the directory specified in the argument, if it doesn't exist.
Returns nothing.
"""
function dir_exists(dir::String)
    try
      readdir(dir)
    catch err
      mkpath(dir)
    end
    return
end

