using BDIO, ALPHAio

@doc raw"""
     open_data_file(path,mode; override=false, comment="")

Open a bdio file in `mode` and return the handle. If a file is created, a path to the file is also made

# Modes
- `"r"`: read mode. Throws an error if the file does not exist
- `"w"`: write mode. if the kwargs `override` is set to `true`, it warns the user and delete any existing file. otherwise throws an error if it exists. `comment`, is passed on to `ALPHAdobs_create`
- `"d"`: delete mode. If open the file and override any existing file without warning the user.  `comment`, is passed on to `ALPHAdobs_create`
 """
function open_data_file(path,mode; comment="",override=false)
    mode == "r" && (return BDIO_open(path,mode))
    if isfile(path)
        if mode == "w"  && override
            @warn "overriding $(basename(path))"
            Base.rm(path)
        elseif mode =="d"
            Base.rm(path,force=true)
        else
            error("$(basename(path)) already exists")
        end
    end
    mkpath(dirname(path))
    return ALPHAdobs_create(path,IOBuffer(comment))
end

@doc """
    close_data_file(fb)

close the data file
"""
close_data_file(fb) = BDIO_close!(fb)

@doc """
     save_data(a,path::String,comment::AbstractString=""; override::Bool = false, extra=nothing)

saves the data `a` in a BDIO file at `path` following the ALPHA conventions. If `override` is given, then
"""
function save_data(a, path::String, comment::AbstractString=""; override::Bool = false, extra=nothing)
    file = open_data_file(path,"w",comment=comment,override=override)
    save_data(a,file,extra=extra)
    ALPHAdobs_close(file)
end

function save_data(a::Dict{String, <:AbstractArray{uwreal}}, fb; extra=nothing)
    kks = collect(keys(a))
    sz =  size(a[kks[1]])
    ALPHAio.ALPHAdobs_vdict_parameters(fb,1, kks, sz; extra=extra)
    for k in values(a)
        ALPHAdobs_add(fb, collect(k))
    end
end

function save_data(a::AbstractArray{<:AbstractArray{uwreal}}, fb; extra=nothing)
    save_data(collect(a[1]),fb,extra=extra)
    for v in a[2:end]
        save_data(collect(v),fb,extra=nothing)
    end
end

save_data(a, fb; extra=nothing) = ALPHAdobs_write(fb,a,extra=extra)

save_data(a;path::String,comment::AbstractString="", override::Bool=false,extra = nothing) = save_data(a,path,comment,override=override, extra=extra)


function read_data(path::AbstractString; get_extra = false)
    file =open_data_file(path,"r")
    res = Any[]
    while ALPHAdobs_next_p(file)
        info =  ALPHAdobs_read_parameters(file)

        type = info["type"]
        nobs = info["nobs"]
        dim  = info["dimensions"]
        size = dim ==0 ?  0 : info["size"]

        if type == "obs"
            if dim == 0
                data = uwreal[]
                for i in 1:nobs
                    push!(data,ALPHAdobs_read_next(file))
                end
            else
                data = [ones(uwreal,size...) for _ in 1:nobs]
                for n in 1:nobs, c in CartesianIndices(data[n])
                    data[n][c] = ALPHAdobs_read_next(file)
                end
            end
        elseif type == "dict"
            keys = info["keys"]
            if dim ==0
                data = [Dict{String,uwreal}(k=>ALPHAdobs_read_next(file)  for k in keys) for _ in 1:nobs]
            else
                data = [Dict{String,Array{uwreal,dim}}() for _ in 1:nobs]
                for n in 1:nobs, k in keys
                    data[n][k] = ones(uwreal,size...)
                    for c in CartesianIndices(data[n][k])
                        data[n][k][c] = ALPHAdobs_read_next(file)
                    end
                end
            end
        end

        if get_extra
            nobs ==1 ? push!(res,(data[1],get(info,"extra",nothing))) : push!(res,(data, get(info,"extra",nothing)))
        else
            nobs ==1 ? push!(res,data[1]) : push!(res,data)
        end
    end

    return length(res) ==1 ?  res[1] : res
end

function save_fit(fit::NamedTuple, path,comment, override=false)
    if !haskey(fit,:par)
        throw(ArgumentError(":par field is missing from fit"))
    end
    _keys = [keys(fit)...] |>x->x[x.!=:par];
    par = getfield(fit,:par)
    extra = Dict(String(k) => getfield(fit,k) for k in _keys)
    save_data(par,path,comment,override=override,extra=extra)
end

function save_fit(fit::T where{T<:AbstractArray{<:NamedTuple}},path,comment,override=false)
    if !all(haskey.(fit,:par))
        throw(ArgumentError(":par field is missing from fit"))
    end
    _keys = keys.(fit) |> x->union([xx for xx in x]...) |> x->filter!(y->y!=:par,x)
    par = getfield.(fit,:par)
    extra = Dict()
    for k in _keys
        extra[String(k)] = [haskey(f,k) ? getfield(f,k) : zeros(Float64) for f in fit]
    end
    save_data(par,path,comment,override=override,extra=extra)
end

function read_fit(path,get_extra = true)
    if !get_extra
        aux = read_data(par)
        return [NamedTuple(:par,aux[i]) for i in eachindex(aux)]
    end
    aux = read_data(path,get_extra=true)
    fit = Vector{NamedTuple}(undef, length(aux))
    extra::Dict = aux[1][2]
    extra["par"] = getfield.(aux,1)
    for i in eachindex(aux)
        fit[i] =NamedTuple((Symbol(key),value[i]) for (key,value) in extra)
    end
    return fit
end

point_to_dict(p::Point) =  Dict(
    "gamma"     => string(p.gamma),
    "x0"        => ismissing(p.x0) ? "moving" : p.x0,
    "qsmearing" => string(p.qsmearing),
    "gsmearing" => string(p.gsmearing)
)

function prop_to_dict(p::Propagator)
    d = Dict()
    d["kappa"] = p.k
    d["mu"] = p.mu
    d["theta"] = [p.theta...]
    d["pF"] = [p.pF...]
    d["src"] = point_to_dict(p.src)
    d["snk"] = point_to_dict(p.snk)
    d["seq_prop"] = p.seq_prop
    return d
end
"""
    function write_corr(C::juobs.Corr; folder=".",ens="ens",set::Union{String,Nothing}=nothing,label::Union{String,Nothing}=nothing,override::Bool=false)

Generates a BDIO file with the correlator `C` using the ALPHA convention (DISCLAIMER: This used the ALPHAio package v0.4.0
written by Alberto Ramos, and as of today the convention is not set to stone, so it may change in future)

The file is saved in `folder` and the naming convention is the following

`ens_set_label_kappa_C.kappa[1]_C.kappa[2]_mu_C.mu[1]_C.mu[2]_tetha1_C.theta1[1]_C.theta1[2]_C.theta1[3]_theta2_C.theta2[1]_C.theta2[2]_C.theta2[3].bdio`

The correlator itself is stored as a vector of uwreal, the correlator parameters (thetas, kappas, mus, y0, gammas) are store in the parameters ui under the `extra` key

if the flag `override` is set to true, it override an existing file with the same name
"""
function write_corr(path,C::Corr;comment = "",override::Bool=false)
    extra = Dict{String,Any}();
    extra["propagators"] = [prop_to_dict(p) for p in C.propagators]
    file = open_data_file(path,"w",override=override,comment = comment)
    ALPHAdobs_write(file,C.obs,extra = extra)
    ALPHAdobs_close(file)
end

function write_corr(C::Corr{2};folder=".",ens="ens",set=nothing,subdirs=nothing,comment="",override::Bool = false )
    dirname = joinpath(folder,ens)
    dirname = isnothing(set) ? dirname : joinpath(dirname,set)
    dirname = joinpath(dirname,"2pt")
    gamma = getfield.(C.points,:gamma) |> x -> join(x,"_")
    dirname = joinpath(dirname,gamma)
    dirname =   isnothing(subdirs) ? dirname : joinpath(dirname,subdirs)
     x0 = getfield.(C.points,:x0) |> x-> join(skipmissing(x),"_")
    kappa = getfield.(C.propagators,:k) |> x-> join(x,"_")
    mu = getfield.(C.propagators,:mu) |> x-> join(x,"_")
    theta1 = join(C.propagators[1].theta,"_")
    theta2 = join(C.propagators[2].theta,"_")

    filename = string(ens,"_x0_",x0,"_",gamma,"_kappa_",kappa,"_mu_",mu,"_theta1_",theta1,
                      "_theta2_",theta2,".bdio")
    write_corr(joinpath(dirname,filename), C,comment=comment, override=override)
end

function write_corr(C::Corr{3};folder=".",ens="ens",set=nothing,subdirs=nothing,comment="",override::Bool = false )
    dirname  = joinpath(folder,ens)
    dirname  = isnothing(set) ? dirname : joinpath(dirname,set)
    dirname  = joinpath(dirname,"3pt")
    gamma    = getfield.(C.points,:gamma) |> x -> join(x,"_")
    dirname  = joinpath(dirname,gamma)
    dirname  = isnothing(subdirs) ? dirname : joinpath(dirname,subdirs)
    x0       = getfield.(C.points,:x0) |> x-> join(skipmissing(x),"_")
    kappa    = getfield.(C.propagators,:k) |> x-> join(x,"_")
    mu       = getfield.(C.propagators,:mu) |> x-> join(x,"_")
    theta1   = join(C.propagators[1].theta,"_")
    theta2   = join(C.propagators[2].theta,"_")
    theta3   = join(C.propagators[3].theta,"_")
    filename = string(ens,"_x0_",x0,"_",gamma,"_kappa_",kappa,"_mu_",mu,"_theta1_",theta1,
                      "_theta2_",theta2,"_theta3_",theta3,".bdio")
    write_corr(joinpath(dirname,filename), C,comment=comment, override=override)
end

dict_to_point(d::Dict) = Point(parse(Gamma,d["gamma"]),
                               d["x0"] == "moving" ? missing : d["x0"],
                               parse(QuarkSmearing.Type,d["qsmearing"]),
                               parse(GluonicSmearing.Type,d["gsmearing"]))

dict_to_prop(d::Dict) = Propagator(d["kappa"],d["mu"],tuple(d["theta"]...),
                                   tuple(d["pF"]...),dict_to_point(d["src"]),
                                   dict_to_point(d["snk"]),d["seq_prop"])
function _read_corr(path)
    obs,extra = read_data(path,get_extra=true)
    prop = tuple((dict_to_prop(d) for d in extra["propagators"])...)
    return Corr(obs,prop)
end

function read_corr(ens;rootdir::String = datadir(),
                   subdir::Union{String,Nothing}=nothing,
                   label::Union{String,Nothing}=nothing,
                   set::Union{String,Nothing}=nothing,
                   theta1::Union{Vector{Float64},Float64} = Float64[],
                   theta2::Union{Vector{Float64},Float64} = Float64[],
                   kappa::AbstractArray{Float64} = Float64[],
                   mu::Vector{Float64} = Float64[])
    dirname = joinpath(rootdir,ens)
    dirname = isnothing(subdir) ? dirname : joinpath(dirname,subdir)
    if !isdir(dirname)
        error("Directory $dirname does not exists")
    end
    files = readdir(dirname)
    !isnothing(set) && filter!(x->contains(x,set),files)
    !isnothing(label) && filter!(x->contains(x,label),files)
    for (spec,txt) in zip([theta1,theta2,kappa,mu],
                          ["theta1","theta2","kappa","mu"])
        length(spec)==0 && continue;
        join([txt;string.(spec)],"_") |>
            aux->filter!(x->contains(x,aux),files)
    end
    if length(files) == 0
        error("No file exist under these conditions")
    end
    if length(files) == 1
        return _read_corr(joinpath(dirname,files[1]))
    else
        return [_read_corr(joinpath(dirname,f)) for f in files]
    end
end
