@enum Gamma Id G0 G1 G2 G3 G5 G0G1 G0G2 G0G3 G0G5 G1G2 G1G3 G1G5 G2G3 G2G5 G3G5   

@doc """
     PointInfo

Store the information of a point in a correlator

# Fields
- `gamma::Gamma`: gamma structure at the point. See also [`Gamma`](@ref)
- `x0::Union{Int64,Missing}`: Time position in lattice units. If `missing` then is the moving point 
"""
struct PointInfo
    gamma::Gamma
    x0::Union{Int64,Missing}
end

struct Propagator
    m::Float64
    mu::Float64
    theta::AbstractVector{Float64}
    pF::AbstractVector{Int64}
    src::Int64; # index to the PointInfo in "insertion"
    snk::Int64; # if 0-> source of the propagator, if -1 sink of the propagator
end

abstract type AbstractCorr end

mutable struct Corr <: AbstractCorr
    obs::AbstractVector
    src::PointInfo
    snk::PointInfo
    insertions::AbstractVector{PointInfo}
    propagators::AbstractVector{Propagator}
    L::Union{Int64,AbstractVector{Int64}} # spatial lenght in lattice units. if only 1 int is given, the lattice volume is L^3, otherwise V = L[1]*L[2]*L[3]
    T::Int64  # Lattice time extent.
    N::Int64 # number of point in the correlator
    function Corr(obs,gammas, xs, ms, mus, thetas,pFs,L)
        N = length(gammas)
        if any(N .!= length.((xs,ms,mus,thetas,pFs,)))
            error("[Corr] Incompatible fields. gammas, xs,ms, mus, thetas and pFs must have the same number of elements")
        end
        src = PointInfo(gammas[1],xs[1])
        snk = PointInfo(gammas[end],xs[end])
        insertions = N == 2 ? PointInfo[] : Vector{PointInfo}(undef,N-2)
        propagators = Vector{Propagator}(undef, N)
        propagators[1] = Propagator(ms[1],mus[1],thetas[1],pFs[1],0, N>2 ? 1 : -1)
        for i in 2:N-1
            xsrc = i-1
            xsnk = i==N-1 ? -1 : i
            propagators[i] = Propagator(ms[i],mus[i],thetas[i],pFs[i],xsrc,xsnk)
            insertions[i-1] = PointInfo(gammas[i+1],xs[i])
        end
        propagators[end] = Propagator(ms[end],mus[end],thetas[end],pFs[end],-1,0)
        T = length(obs)
        return new(obs,src,snk,insertions,propagators,L,T,N)
    end   
end

import Base:show

function show_customized_PointInfo(io::IO, p::PointInfo;tab="")
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
    println(io,tab,"(m,mu): (",p.m," ",p.mu,")")
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
    show_customized_PointInfo(io,c.src,tab="   ");
    if (N!=2)
        for (i,p) in pairs(c.insertions)
            println("$(i) insertion point:")
            show_customized_PointInfo(io,p,tab="    ")
        end
    end
    println("sink:")
    show_customized_PointInfo(io,c.snk,tab="    ")
    label = ["$(i) insersition point" for i in 1:N-2]
    for (i,p) in pairs(c.propagators)
        println("propagator $(i):")
        show_customized_propagator(io,p,label,tab="    ")
    end
end


