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
    Point() = new(None,-1,None,None)
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
    src::Int64
    snk::Int64
    Propagator() = new(0.,0.,(0.,0.,0),(0,0,0,0),-1,-1)
    Propagator(k,m,t,p,src,snk) = new(k,m,t,p,src,snk)
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
    Corr(o,po,pr) = new{N}(o,po,pr)
    Corr() = new{N}([],ntuple(x->Point(),N),ntuple(x->Propagator(),N))
    Corr(x::Base.Generator) = new{N}(x...)
end

function update(c::Corr{N} where N; k...)
    isempty(k) && (return c)
    return Corr(get(k,s,getfield(c,s)) for s in fieldnames(Corr{N}))
end

import Base:show

function show_customized_Point(io::IO, p::Point;tab="")
    println(io,tab,"gamma: ",p.gamma)
    println(io,tab,"x0:    ",ismissing(p.x0) ? "varying" : p.x0)
end

function show_customized_propagator(io::IO, p::Propagator,label::Vector{<:AbstractString};tab="")
    function f(x)
        if x ==0
            return "source"
        end
        if x==-1
            return "sink"
        else
            return label[x]
        end
    end
    println(io,tab,"(k,mu): (",p.k," ",p.mu,")")
    println(io,tab,"pF:     [",join(string.(p.pF),", "),"]")
    println(io,tab,"theta:  [",join(string.(p.theta),", "),"]")
    println(io,tab,"propagate from ", f(p.src)," to ",f(p.snk))
end

function show(io::IO,c::Corr)
    println(io,"$(c.N)-point correlator:")
    if c.L isa Int64
        println(io,"Lattice size = $(c.L)^3*$(c.T)")
    else
        println(io,"Lattice size = $(c.L[1])*$(c.L[2])*$(c.L[3])*$(c.T)")
    end
    println("source: ")
    show_customized_Point(io,c.src,tab="   ");
    if (N!=2)
        for (i,p) in pairs(c.insertions)
            println("$(i) insertion point:")
            show_customized_Point(io,p,tab="    ")
        end
    end
    println("sink:")
    show_customized_Point(io,c.snk,tab="    ")
    label = ["$(i) insersition point" for i in 1:N-2]
    for (i,p) in pairs(c.propagators)
        println("propagator $(i):")
        show_customized_propagator(io,p,label,tab="    ")
    end
end

function show(io::IO,p::Point)
    print(io,"gamma = ", p.gamma, ",\tx0 = ", p.x0)
    print(io,",\tquark smering = ",p.qsmearing,",\tgluonic smearing = ",p.gsmearing)
end

show(io::IO,p::Propagator) = show_customized_propagator(io,p,["",""])

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
