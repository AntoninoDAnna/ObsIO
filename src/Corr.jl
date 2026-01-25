import Base:show

@doc raw"""
    struct Point

Immutable structure that contains all information regarding a point in a correlator

## Fields
- `gamma::Gamma`: gamma structure (See also [`Gamma`](@ref))
- `x0::Union{Int64,Missing}`: Position in lattice unit. If `missing` then it rappresent the varying point.
- `qsmearing::QuarkSmearing.Type`: Quark Smearing at x0 (See also [`QuarkSmearing](@ref))
- `gsmearing::GluonicSmearing.Type`: Gluonic Smearing at x0 (See also [`GluonicSmearing](@ref))
        """
struct Point
    gamma::Gamma
    x0::Union{Int64,Missing}
    qsmearing::QuarkSmearing.Type
    gsmearing::GluonicSmearing.Type
    Point() = new(None,-1,QuarkSmearing.None,GluonicSmearing.None)
    Point(g,x,qs,gs) = new(g,x,qs,gs)
    Point(x::Base.Generator) =new(x...)
end

@doc raw"""
    struct Propagator

Immutable structure that contains all the informations regarding a propagator

## Fields

- `k::Float64, mu::Float64`: hopping parameter and twisted mass of the quark
- `theta::NTuple{3,Float64}`: Twisted Boundary Condition associated to the propagator that induce momentum on the quark
- `pF::NTuple{4,Int64}`: Quantized momentum associate to the propagator.
- `src::Point,snk::Point`: source and sink Point (See also [`Point`](@ref))
- `seq_prop::Bool`: flag for sequential propagator.

## Warning
  `seq_prop` accept also `Int64`. This is an internal functionality used to keep information on the source propagator and should not be used externally.
    """
struct Propagator
    k::Float64
    mu::Float64
    theta::NTuple{3,Float64}
    pF::NTuple{4,Int64}
    src::Point
    snk::Point
    seq_prop::Union{Int64,Bool}
    Propagator() = new(0.,0.,(0.,0.,0),(0,0,0,0),Point(),Point(),-1)
    Propagator(k,m,t,p,src,snk,seq_prop=false) = new(k,m,t,p,src,snk,seq_prop)
    Propagator(x::Base.Generator) = new(x...)
end

@doc raw"""
     AbstractCorr

Abstract Type for correlators
"""
abstract type AbstractCorr end

@doc raw"""
     struct Corr{N} <:AbstractCorr

immutable structure for `N` point correlators.
This structure support correlators where only one `Point` is
moving and the other `N-1` point are fixed.

## Fields

- `obs::AbstractVector`: Observable Vector. It contains the correlator data. Time dependande is coded in the index.
- `points::NTuple{N,Point}`: Correlator's Point. By convention, `Point`s are ordered in sequence with the first point being the source and the last point the sink (See also [`Point`](@ref))
- `propagators::NTuple{N,Proopagator`: Propagators list. By convention, `Propagator`s are ordered in sequence. The sink of `propagator[i]` is the source of `propagators[i+1]`. Moreover, the source of `propagator[i]` is `points[i]` and the sink is `point[i%N +1]`. (See also [`Propagator`](@ref))
"""
struct Corr{N} <: AbstractCorr
    obs::AbstractVector
    points::NTuple{N,Point}
    propagators::NTuple{N,Propagator}
    function Corr(o,pr::NTuple{N,Propagator}) where N
        pts = (pr[1].src,)
        for i in 1:N-1
            pts = (pts..., pr[i].snk)
        end
        return new{N}(o,pts,pr)
    end
    Corr(o,po::NTuple{N,Point},pr::NTuple{N,Propagator}) where N=
        new{N}(o,po,pr)
    Corr(N::Int64) = new{N}([],ntuple(x->Point(),N),ntuple(x->Propagator(),N))
    Corr(x::Base.Generator) = Corr(x...)
end

@doc raw"""
     kappa(c::Corr{N}) where N

return a tuple containing the hopping parameters of the propagator
"""
kappa(c::Corr{N}) where N = ntuple(x->c.propagator[x].k, N)

@doc raw"""
     mu(c::Corr{N}) where N

return a tuple containing the twisted masses of the propagator
"""
mu(c::Corr{N}) where N = ntuple(x->c.propagator[x].mu, N)

@doc raw"""
     theta(c::Corr{N}) where N

return a tuple containing the theta boundary conditions of the propagator
"""
theta(c::Corr{N}) where N = ntuple(x->c.propagator[x].theta, N)

@doc raw"""
     src(c::Corr{N}) where N

return the source position in lattice units. Equivalent to `c.point[1].x0`
"""
src(c::Corr{N}) where N = c.point[1].x0

@doc raw"""
     src(c::Corr{N}) where N

return the source position in lattice units. Equivalent to `c.point[end].x0`
"""
snk(c::Corr{N}) where N = c.point[end].x0

@doc raw"""
     src(c::Corr{N}) where N

return the source-sink distance in lattice. assume that `N ⪩ 3`
    """
ts(c::corr{N}) where N = snk(c) - src(c)


@doc raw"""
     __update__(p::Point;k...)
     __update__(p::Propagator;k...)
     __update__(c::Corr{N} where N; k...)

return a new object that updates the original object with the values in `k`. If no new information is given, then return the original

## Warning

This function is an internal method used to simplify the reading data files process. It is not meant to be used by the external users.
"""
function __update__(p::Point;k...)
    isempty(k) && (return p)
    return Point(get(k,s,getfield(p,s)) for s in fieldnames(Point))
end

function __update__(p::Propagator;k...)
    isempty(k) && (return p)
    return Propagator(get(k,s,getfield(p,s)) for s in fieldnames(Propagator))
end

function __update__(c::Corr{N} where N; k...)
    isempty(k) && (return c)
    return Corr(get(k,s,getfield(c,s)) for s in fieldnames(Corr))
end

@doc raw"""
     struct GlobalHeader

Contains all the information inside the global header of a meson.dat file.
"""
struct GlobalHeader
    ncorr::Int32
    nnoise::Int32
    tvals::Int32
    noise::Noise.Type
    hsize::Int32
    GlobalHeader(a,b,c,d) = new(a,b,c,d,4*4)
    GlobalHeader(x::Vector{Int32}) = new(x[1],x[2],x[3],Noise.Type(x[4]),4*4)
end

@doc raw"""
     mutable struct Smearing{T<:EnumClass}

contains the Smearing informations. See also [`GluonicSmearing`](@ref), [`QuarkSmearing`](@ref)
"""
mutable struct Smearing{T<:EnumClass}
    type::T
    niter::Int32
    eps::Float64
    Smearing(t::T) where {T<:EnumClass} = new{T}(t, 0, 0.0)
    Smearing(t::T, n, e,) where {T<:EnumClass} = new{T}(t, n, e)
end

@doc raw """
     mutable struct CorrHeader

contains the information inside a correlator Header of a meson.dat file
"""
mutable struct CorrHeader
    k::NTuple{2,Float64}
    mu::NTuple{2,Float64}
    dp::NTuple{2,Float64}
    theta::Vector{Vector{Float64}}
    q::NTuple{2,Smearing}
    g::NTuple{2,Smearing}
    type::NTuple{2,Gamma}
    x0::Int32
    is_real::Bool
    hsize::Int32 #headersize
    dsize::Int32 #datasize / (nnoise * T * ncfg)
    function CorrHeader(aux_f::Vector{Float64}, aux_i::Vector{Int32}, theta::Vector{Float64}, sm_par::Vector{Smearing})
        a = new()
        a.k  = (aux_f[1],aux_f[2])
        a.mu = (aux_f[3],aux_f[4])
        a.dp = (aux_f[5],aux_f[6])
        a.type = Gamma.((aux_i[1],aux_i[2]))
        a.x0 = aux_i[3]
        a.is_real = aux_i[4]==1
        a.theta = [theta[1:3],theta[4:6]]
        a.q = (sm_par[1],sm_par[2])
        a.g = (sm_par[3],sm_par[4])
        a.hsize = 8*12 + 4*8 #without smearing parameters niter, neps
        for q in a.q
            q.type == QuarkSmearing.Local && continue;
            a.hsize += 8+4
        end
        for g in a.g
            if !(g.type == GluonicSmearing.APE || g.type == GluonicSmearing.WilsonFlow3D)
                continue;
            end
            a.hsize += 8+4
        end
        a.dsize = a.is_real ? 8 : 16
        return a
    end
    function CorrHeader(aux_f::Vector{Float64}, aux_i::Vector{Int32})
        a = new()
        a.k  = (aux_f[1],aux_f[2])
        a.mu = (aux_f[3],aux_f[4])
        a.dp = zeros(2)
        a.type = Gamma.((aux_i[1],aux_i[2]))
        a.x0 = aux_i[3]
        a.is_real = aux_i[4]==1
        a.theta= [zeros(3),zeros(3)]
        a.hsize = 8*12 + 4*8 #without smearing parameters niter, neps
        a.dsize = 16 - 8* a.is_real
        return a
    end
end

@doc raw"""
     mutable struct CorrData

Contains all the correlator information saved inside a meson.dat
"""
mutable struct CorrData
    header::CorrHeader
    vcfg::Array{Int32}
    re_data::Array{Float64}
    im_data::Array{Float64}
    id::String
    CorrData(a, b, c, d, e) = new(a, b, c, d, e,)
    CorrData(header, ncfg, nt, id) = new(header,collect(1:ncfg),
                                         zeros(ncfg,nt),zeros(ncfg,nt),id)
end

function Base.show(io::IO, s::Smearing)
    ST = typeof(s.type)
    if isa(s.type,QuarkSmearing.Type) && s.type == QuarkSmearing.Local
        print(io,"Quark Smering: ",s.type)
    elseif isa(s.type,QuarkSmearing.Type)
        print(io,"Quark Smering: ",s.type,", niter: ",s.niter,", eps: ",s.eps)
    elseif isa(s.type,GluonicSmearing.Type) && s.type == GluonicSmearing.Local
        print(io,"Gluonic Smering: ",s.type)
    elseif isa(s.type,GluonicSmearing.Type)
        print(io,"Gluonic Smering: ",s.type,", niter: ",s.niter,", eps: ",s.eps)
    end
end

function Base.show(io::IO, c::CorrHeader)
    println(io,"kappa: ",c.k,". ")
    println(io,"mu: ",c.mu,". ")
    println(io,"dp: ", c.dp,". ")
    println(io,"theta: ",c.theta[1],", ",c.theta[2],". ")
    println(io,c.q[1],", ",c.q[2],", ")
    println(io,c.g[1],", ",c.g[2],", ")
    println(io,"type: ",c.type,", ")
    print(io,"source: ", c.x0,", ")
end

function show(io::IO, c::CorrData)
    println(io, "Correlator data. Ensemble ",c.id)
    println(io, "Configuration length ", c.vcfg)
    println(io, c.header)
end

function show_customized_propagator(io::IO, p::Propagator,label::Vector{<:AbstractString};tab="")
    print(io,tab,"(k,mu): (",p.k," ",p.mu,")")
    print(io,tab,"pF:     [",join(string.(p.pF),", "),"]")
    print(io,tab,"theta:  [",join(string.(p.theta),", "),"]")
    print(io,tab,"propagate from ", f(p.src)," to ",f(p.snk))
end

function show(io::IO,c::Corr{N} where N)
    N = length(c.points)
    println(io, "$(N)-point correlator")
    println(io,"Points:")
    for p in c.points
        println(io,"\t",p,)
        println(io,"-"^80)
    end
    println(io,"\n","Propagators:")
    for p in c.propagators
        println(io, "\t", p)
        println(io,"-"^80)
    end
end

function show(io::IO,p::Point)
    print(io,"gamma = ", p.gamma, ",\tx0 = ", p.x0)
    print(io,",\tquark smering = ",p.qsmearing,",\tgluonic smearing = ",p.gsmearing)
end

function show(io::IO,p::Propagator)
    print(io, "kappa = ", p.k," ")
    print(io, "mu = ", p.mu, " ")
    print(io, "pF = ", p.pF, " ")
    println(io, "theta = ",p.theta," ")
    if isa(p.seq_prop, Bool) && p.seq_prop
        println(io, "It's a sequential propagator")
    elseif isa(p.seq_prop,Int64)
        println(io, "It's a sequential propagator starting from propagator ",
                p.seq_prop)
    end
    println(io, "source: \n", p.src)
    println(io, "sink: \n", p.snk)
end
