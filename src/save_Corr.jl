point_to_dict(p::Point{<:Integer}) =  Dict(
    "gamma"     => string(p.gamma),
    "x0"        => string(p.x0),
    "qsmearing" => string(p.qsmearing),
    "gsmearing" => string(p.gsmearing)
)

point_to_dict(p::Point{Missing}) =  Dict(
    "gamma"     => string(p.gamma),
    "x0"        => "moving",
    "qsmearing" => string(p.qsmearing),
    "gsmearing" => string(p.gsmearing)
)

str2bc = Dict(split(string(T),".")[end]=>T for T in (Open,SF,Open_SF,Periodic))


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
    function write_corr(C::Corr; folder=".",ens="ens",set::Union{String,Nothing}=nothing,label::Union{String,Nothing}=nothing,override::Bool=false)

Generates a BDIO file with the correlator `C` using the ALPHA convention (DISCLAIMER: This used the ALPHAio package v0.4.0
written by Alberto Ramos, and as of today the convention is not set to stone, so it may change in future)

The file is saved in `folder` and the naming convention is the following

`ens_set_label_kappa_C.kappa[1]_C.kappa[2]_mu_C.mu[1]_C.mu[2]_tetha1_C.theta1[1]_C.theta1[2]_C.theta1[3]_theta2_C.theta2[1]_C.theta2[2]_C.theta2[3].bdio`

The correlator itself is stored as a vector of uwreal, the correlator parameters (thetas, kappas, mus, y0, gammas) are store in the parameters ui under the `extra` key

if the flag `override` is set to true, it override an existing file with the same name
"""
function write_corr(path,C::Corr{N,BC,T};comment = "",override::Bool=false) where {N, BC<:AbstractBC,T}
    extra = Dict{String,Any}();
    extra["propagators"] = [prop_to_dict(p) for p in C.propagators]
    extra["boundary conditions"] = split(string(BC),".")[end]
    file = open_data_file(path,"w",override=override,comment = comment)
    ALPHAdobs_write(file,C.obs,extra = extra)
    ALPHAdobs_close(file)
end

function write_corr(C::Corr{2,BC,T}; folder=".",ens="ens",set=nothing,subdirs=nothing,comment="",override::Bool = false ,info=false) where {BC<:AbstractBC,T}
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
    info && println(joinpath(dirname,filename))
    write_corr(joinpath(dirname,filename), C,comment=comment, override=override)
end

function write_corr(C::Corr{3,BC,T};folder=".",ens="ens",set=nothing,subdirs=nothing,comment="",override::Bool = false, info=false ) where {BC <:AbstractBC, T}
    dirname  = joinpath(folder,ens)
    dirname  = isnothing(set) ? dirname : joinpath(dirname,set)
    dirname  = joinpath(dirname,"3pt")
    gamma    = getfield.(C.points,:gamma) |> x -> join(x,"_")
    dirname  = joinpath(dirname,gamma)
    dirname  = isnothing(subdirs) ? dirname : joinpath(dirname,subdirs)
    x0       = getfield.(C.points,:x0) |> x-> ºjoin(skipmissing(x),"_")
    kappa    = getfield.(C.propagators,:k) |> x-> join(x,"_")
    mu       = getfield.(C.propagators,:mu) |> x-> join(x,"_")
    theta1   = join(C.propagators[1].theta,"_")
    theta2   = join(C.propagators[2].theta,"_")
    theta3   = join(C.propagators[3].theta,"_")
    filename = string(ens,"_x0_",x0,"_",gamma,"_kappa_",kappa,"_mu_",mu,"_theta1_",theta1,
                      "_theta2_",theta2,"_theta3_",theta3,".bdio")
    info &&println(joinpath(dirname,filename))
    write_corr(joinpath(dirname,filename), C,comment=comment, override=override)
end
