@doc """
     Point

Store the information of a point in a correlator

# Fields
- `gamma::Gamma`: gamma structure at the point. See also [`Gamma`](@ref)
- `x0::Union{Int64,Missing}`: Time position in lattice units. If `missing` then is the moving point
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

function update(p::Point;k...)
    isempty(k) && (return p)
    return Point(get(k,s,getfield(p,s)) for s in fieldnames(Point))
end

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

function update(p::Propagator;k...)
    isempty(k) && (return p)
    return Propagator(get(k,s,getfield(p,s)) for s in fieldnames(Propagator))
end

abstract type AbstractCorr end

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

function update(c::Corr{N} where N; k...)
    isempty(k) && (return c)
    return Corr(get(k,s,getfield(c,s)) for s in fieldnames(Corr))
end

import Base:show

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

@doc"""
        read_data(path;kwargs...)

It read the data in `path` according to `kwargs`
This function call the relevant function according to the extention of `path`

"""
function read_data(path; kwargs...)
    if !isfile(path)
        error("file $path does not exist")
    end
    b = basename(path)
    if contais("mesons.dat", "")  == join(b[end-1:end],".")
        return read_mesons(path;kwargs...)
    else
        error("This type of file is not supported")
    end
end
