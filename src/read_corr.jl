struct NoMatchingFileError <:Exception
    message::String
end

function Base.showerror(io::IO,e::NoMatchingFileError)
    print("NoMatchingFileError: ", e.message)
end

to_int(x::Int64) = x
to_int(x::AbstractString) = parse(Int64,x)
dict_to_point(d::Dict) = ObsIO.Point(parse(ObsIO.Gamma,d["gamma"]),
                               d["x0"] == "moving" ? missing : to_int(d["x0"]),
                               parse(ObsIO.QuarkSmearing.Type,d["qsmearing"]),
                               parse(ObsIO.GluonicSmearing.Type,d["gsmearing"]))

dict_to_prop(d::Dict) = ObsIO.Propagator(d["kappa"],d["mu"],tuple(d["theta"]...),
                                   tuple(d["pF"]...),dict_to_point(d["src"]),
                                   dict_to_point(d["snk"]),d["seq_prop"])


function read_bc(d::Dict)
    s = d["boundary conditions"]
    return get(str2bc,s,Open)
end


function _read_corr(path)
    obs,extra = read_data(path,get_extra=true)
    prop = tuple((dict_to_prop(d) for d in extra["propagators"])...)
    bc = read_bc(extra)
    return Corr(obs,prop,bc)
end



function make_filter(;filters...)
    filter(x::String) = all(contains(x,v) for (_,v) in filters)
    return filter
end


nofile_found(dirname;filters...) =  throw(NoMatchingFileError(string("No file found in ", dirname, " that fulfills the requirements: ", join(["$k => $v" for (k,v) in filters], "; "))))

empty_folder(dirname) = throw(NoMatchingFileError(string("Folder ", dirname, " is empty")))

function __find_corr_file(ens::String; rootdir::String,
                          subdir::String = "",
                          filters...)
    dirname = joinpath(rootdir,ens,subdir)
    isdir(dirname) || error("$dirname does not exists")
    files = readdir(dirname)
    isempty(files) && empty_folder(dirname)
    files = filter(make_filter(;filters...),files)
    isempty(files) && nofile_found(dirname;filters...)
    return joinpath.(dirname,files)
end

@doc raw"""
    read_corr(ens::String; rootdir::String = ".", subdir::String = "", filters...)

Read all the correlators in folder `joinpath(rootdir,ens,subdir)`. If `filters` are provided, it read the correlators which respect the filters. The filters act on the filenames as `contains(filename,filter)`, so it can be any object accepted by the function `contains` as a second arguments. If no file in found, an error is thrown

"""
function read_corr(ens;rootdir::String = ".", subdir::String = "", filters...)
    files = __find_corr_file(ens,rootdir=rootdir, subdir=subdir; filters...)
    if length(files) == 1
        return _read_corr(files[1])
    else
        return [_read_corr(f) for f in files]
    end
end
