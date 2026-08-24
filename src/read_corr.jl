dict_to_point(d::Dict) = ObsIO.Point(parse(ObsIO.Gamma,d["gamma"]),
                               d["x0"] == "moving" ? missing : d["x0"],
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


function nofile_found(dirname;filters...)
   error(LazyString(" No file found in $dirname that fulfills the requirements: ", join(["$k => $v" for (k,v) in filters], "; ")))
end

function __find_corr_file(ens::String; rootdir::String,
                          subdir::String = "",
                          filters...)
    dirname = joinpath(rootdir,ens,subdir)
    isdir(dirname) || error("$dirname does not exists")
    files = filter(make_filter(;filters...),readdir(dirname,join=true))
    !isempty(files) || nofile_found(dirname;filters...)
    return files
end


function read_corr(ens;rootdir::String = datadir(),
                   subdir::String = "",
                   filters...)
    files = __find_corr_file(ens,rootdir=rootdir, subdir=subdir; filters...)
    if length(files) == 1
        return _read_corr(files[1])
    else
        return [_read_corr(f) for f in files]
    end
end
