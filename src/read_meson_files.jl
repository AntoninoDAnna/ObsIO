using Statistics, ADerrors
struct GlobalHeader
    ncorr::Int32
    nnoise::Int32
    tvals::Int32
    noise::Noise.Type
    hsize::Int32
    GlobalHeader(a,b,c,d) = new(a,b,c,d,4*4)
    GlobalHeader(x::Vector{Int32}) = new(x[1],x[2],x[3],Noise.Type(x[4]),4*4)
end

function read_GlobalHeader(path::String)
    data = open(path, "r")
    g_header = zeros(Int32, 4)
    read!(data, g_header)
    close(data)
    a = GlobalHeader(g_header)
    return a
end

mutable struct Smearing{T<:EnumClass}
    type::T
    niter::Int32
    eps::Float64
    Smearing(t::T) where {T<:EnumClass} = new{T}(t, 0, 0.0)
    Smearing(t::T, n, e,) where {T<:EnumClass} = new{T}(t, n, e)
end

import Base: show
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

function read_CorrHeader(path::String; legacy::Bool=false)
    gh = read_GlobalHeader(path)
    data = open(path, "r")
    seek(data, gh.hsize)
    a = Vector{CorrHeader}(undef, gh.ncorr)
    if !legacy
        aux_f = zeros(Float64, 6)
        aux_i = zeros(Int32, 4)
        theta = zeros(Float64, 6)
        for k = 1:gh.ncorr
            read!(data, aux_f)
            read!(data, theta)
            qs1 = read(data, Int32) |> QuarkSmearing.Type
            if qs1 != QuarkSmearing.Local
                qn1 = read(data, Int32)
                qeps1 = read(data, Float64)
                q1 = Smearing(qs1, qn1, qeps1)
            else
                q1 = Smearing(qs1)
            end
            qs2 = read(data, Int32)  |> QuarkSmearing.Type
            if qs2 != QuarkSmearing.Local
                qn2 = read(data, Int32)
                qeps2 = read(data, Float64)
                q2 = Smearing(qs2, qn2, qeps2)
            else
                q2 = Smearing(qs2)
            end
            gs1 = read(data, Int32)
            if gs1 != 0 && gs1 != 3 && gs1 != 4
                gn1 = read(data, Int32)
                geps1 = read(data, Float64)
                g1 = Smearing(GluonicSmearing.Type(gs1), gn1, geps1)
            elseif gs1 == 3 || gs1 == 4
                g1 = Smearing(GluonicSmearing.Type(gs1), q1.niter, q1.eps)
            else
                g1 = Smearing(GluonicSmearing.Type(gs1))
            end
            gs2 = read(data, Int32)
            if gs2 != 0 && gs2 != 3 && gs2 != 4
                gn2 = read(data, Int32)
                geps2 = read(data, Float64)
                g2 = Smearing(GluonicSmearing.Type(gs2), gn2, geps2)
            elseif gs1 == 3 || gs1 == 4
                g2 = Smearing(GluonicSmearing.Type(gs2), q2.niter, q2.eps)
            else
                g2 = Smearing(GluonicSmearing.Type(gs2))
            end
            read!(data, aux_i)
            a[k] = CorrHeader(aux_f, aux_i, theta, [q1, q2, g1, g2])
        end
    else
        aux_f = zeros(Float64, 4)
        aux_i = zeros(Int32, 4)
        for k = 1:gh.ncorr
            read!(data, aux_f)
            read!(data, aux_i)
            a[k] = CorrHeader(aux_f, aux_i)
        end
    end
    close(data)
    return a
end

mutable struct CorrData
    header::CorrHeader
    vcfg::Array{Int32}
    re_data::Array{Float64}
    im_data::Array{Float64}
    id::String
    CorrData(a, b, c, d, e) = new(a, b, c, d, e,)
end

function show(io::IO, c::CorrData)
    println(io, "Correlator data. Ensemble ",c.id)
    println(io, "Configuration length ", c.vcfg)
    println(io, c.header)
end

function find_match(ch,g1::Gamma,g2::Gamma)
    function f(x)
        return (x.type[1] == g1 || g1 ==None) &&  (x.type[2] == g2 || g2==None)
    end
    findall(f,ch)
end

function find_match(ch,gamma::NTuple{2,Gamma}...)
    function f(x)
        return x.type in gamma
    end
    findall(f,ch)
end

function _read_mesons(path,gh,ch,match;nnoise_trunc,legacy,id, correction::Bool = false)
    ncorr = gh.ncorr
    tvals = gh.tvals
    nnoise = gh.nnoise
    nnoise_trunc = isnothing(nnoise_trunc) ? nnoise : min(nnoise, nnoise_trunc)
    ncfg = let
        fsize = filesize(path)
        datsize = 4 + sum(getfield.(ch, :dsize)) * tvals * nnoise #data_size / ncnfg
        div(fsize - gh.hsize - sum(getfield.(ch, :hsize)), datsize)
        #(total size - header_size) / data_size
    end
    N = correction ? div(length(match),2) : length(match)
    data_re = zeros(N, ncfg, tvals)
    data_im = zeros(N, ncfg, tvals)
    vcfg = Array{Int32}(undef, ncfg)
    open(path,"r") do data
        seek(data, gh.hsize + sum(c.hsize for c in ch))
        for icfg = 1:ncfg
            vcfg[icfg] = read(data, Int32)
            c,sgn=1,+1 ## if correction, at div(ncorr,2) it changes to -
            for k = 1:ncorr
                if !(k in match)
                    seek(data, position(data)  + ch[k].dsize*tvals*nnoise)
                    continue
                end
                if ch[k].is_real
                    tmp = Array{Float64}(undef, tvals*nnoise)
                    read!(data, tmp)
                    tmp2 = reshape(tmp, (nnoise, tvals))
                    tmp2 = mean(tmp2[1:nnoise_trunc, :], dims=1)
                    data_re[c, icfg, :] .+= sgn * tmp2[1, :]
                else
                    tmp = Array{Float64}(undef, 2*tvals*nnoise)
                    read!(data, tmp)
                    tmp2 = reshape(tmp, (2, nnoise, tvals))
                    tmp2 = mean(tmp2[:, 1:nnoise_trunc, :], dims=2)
                    data_re[c, icfg, :] .+= sgn*tmp2[1, 1, :]
                    data_im[c, icfg, :] .+= sgn*tmp2[2, 1, :]
                end
                if k==div(ncorr,2)
                    c,sgn = 1,-1
                else
                    c+=1
                end
            end
        end
    end
    res = Array{CorrData}(undef, N)
    for (c,cm) in pairs(match)
        res[c] = CorrData(ch[cm], vcfg, data_re[c, :, :], data_im[c, :, :], id)
    end
    return res
end

function get_id(path,id)
    !isnothing(id) && return id
    bname = basename(path)
    m = findfirst(r"[A-Z][0-9]{3}r[0-9]{3}", bname)
    return bname[m[1:4]]
end

@doc raw"""

    read_mesons(path::String, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false)

    read_mesons(path::Vector{String}, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false)

This function read a mesons dat file at a given path and returns a vector of `CData` structures for different masses and Dirac structures.
Dirac structures `g1` (source) and/or `g2` (sink) can be passed as string arguments in order to filter correaltors.
ADerrors id can be specified as argument. If is not specified, the `id` is fixed according to the ensemble name (example: "H400"-> id = "H400")

*For the old version (without smearing, distance preconditioning and theta) set legacy=true.

Examples:
```@example
read_mesons(path)
read_mesons(path, "G5")
read_mesons(path, nothing, "G5")
read_mesons(path, "G5", "G5")
read_mesons(path, "G5", "G5", id="H100")

```
    """
function read_mesons(path::String,
                     g1::Gamma = None,
                     g2::Gamma = None;
                     id::Union{String, Nothing}=nothing,
                     legacy::Bool=false,
                     nnoise_trunc::Union{Int64, Nothing}=nothing)
    id = get_id(path,id)
    g_header = read_GlobalHeader(path)
    c_header = read_CorrHeader(path, legacy=legacy)
    corr_match = find_match(c_header,g1,g2)
    return  _read_mesons(path,g_header,c_header,corr_match,
                         nnoise_trunc = nnoise_trunc, legacy=legacy,id=id)
end

function read_mesons(path::String,
                     gamma::NTuple{2,Gamma}...;
                     id::Union{String, Nothing}=nothing,
                     legacy::Bool=false,
                     nnoise_trunc::Union{Int64, Nothing}=nothing)
    id = get_id(path,id)
    g_header = read_GlobalHeader(path)
    c_header = read_CorrHeader(path, legacy=legacy)
    corr_match = find_match(c_header,gamma...)
    return _read_mesons(path,g_header,c_header,corr_match,
                         nnoise_trunc = nnoise_trunc, legacy=legacy,id=id)
end

function read_mesons(path::Vector{String}, p...;k... )
    res = [read_mesons(_path,p...;k...) for _path in path]
    nrep = length(res)
    ncorr = length(res[1])
    cdata = Vector{Vector{CorrData}}(undef, ncorr)
    for icorr = 1:ncorr
        cdata[icorr] = Vector{CorrData}(undef, nrep)
        for r = 1:nrep
            cdata[icorr][r] = res[r][icorr]
        end
    end
    return cdata
end

function read_mesons_correction(path::String,
                                g1::Gamma = None,
                                g2::Gamma = None;
                                id::Union{String, Nothing}=nothing,
                                legacy::Bool=false,
                                nnoise_trunc::Union{Int64, Nothing}=nothing)
    id = get_id(path,id)
    g_header = read_GHeader(path)
    c_header = read_CHeader(path, legacy=legacy)
    match = find_match(c_header,g1,g2)
    return _read_mesons(path,g_header,c_header,match,id=id,nnoise_trunc=nnoise_trunc,
                        legacy=legacy,correction=true)
end

function read_mesons_correction(path::String,
                                gamma::NTuple{2,Gamma}...;
                                id::Union{String, Nothing}=nothing,
                                legacy::Bool=false,
                                nnoise_trunc::Union{Int64, Nothing}=nothing)
    id = get_id(path,id)
    g_header = read_GHeader(path)
    c_header = read_CHeader(path, legacy=legacy)
    match = find_match(c_header,gamma...)
    return _read_mesons(path,g_header,c_header,match,id=id,nnoise_trunc=nnoise_trunc,
                        legacy=legacy,correction=true)
end

function read_mesons_correction(path::Vector{String}, p...;k... )
    res = [read_mesons_correction(_path,p...;k...) for _path in path]
    nrep = length(res)
    ncorr = length(res[1])
    cdata = Vector{Vector{CorrData}}(undef, ncorr)
    for icorr = 1:ncorr
        cdata[icorr] = Vector{CData}(undef, nrep)
        for r = 1:nrep
            cdata[icorr][r] = res[r][icorr]
        end
    end
    return cdata
end

function apply_rw(data::Array{Float64}, W::Matrix{Float64}, vcfg::Union{Nothing, Vector{Int32}}=nothing; id::Union{String, Nothing}=nothing, fs::Bool=false)
    nc = size(W,2)
    if isnothing(vcfg)
        vcfg = collect(1:nc)
    end
    if fs == false
        if id == "A653"
            rw1 = W[1, 1:nc]
            rw = rw1
            data_r = data .* rw[vcfg]
            return (data_r, rw)
        else
            rw1 = W[1, 1:nc]
            rw2 = W[2, 1:nc]
            rw = rw1 .* rw2
            data_r = data .* rw[vcfg]
            return (data_r, rw)
        end
    else
        rw1 = W[1, 1:nc]
        rw2 = W[2, 1:nc]
        rw_s = [1.0 for i in 1:nc]
        if id in ["H105r001", "H105r002", "H105r005", "J303", "J303r003"]
            rw_s[flag_s[id]] .= -1.0
        end
        rw = rw1 .* rw2 .* rw_s
        data_r = data .* rw[vcfg]
        return (data_r, rw)
    end
end

function apply_rw(data::Vector{<:Array{Float64}}, W::Vector{Matrix{Float64}}, vcfg::Union{Nothing, Vector{Vector{Int32}}}=nothing; id::Union{String, Nothing}=nothing, fs::Bool=false)
    nc = size.(W, 2)
    if isnothing(vcfg)
        vcfg = [collect(1:nc[k]) for k=1:length(nc)]
    end
    if fs == false
        if id == "A653"
            rw1 = [W[k][1, 1:nc[k]] for k=1:length(W)]
            rw = [rw1[k] for k =1:length(W)]
            data_r = [data[k] .* rw[k][vcfg[k]] for k=1:length(data)]
            return (data_r, rw)
        else
            rw1 = [W[k][1, 1:nc[k]] for k=1:length(W)]
            rw2 = [W[k][2, 1:nc[k]] for k=1:length(W)]
            rw = [rw1[k] .* rw2[k] for k =1:length(W)]
            data_r = [data[k] .* rw[k][vcfg[k]] for k=1:length(data)]
            return (data_r, rw)
        end
    else
        rw1 = [W[k][1, 1:nc[k]] for k=1:length(W)]
        rw2 = [W[k][2, 1:nc[k]] for k=1:length(W)]
        rw_s = [[1.0 for i in 1:nc[k]] for k=1:length(W)]
        if id == "H105"
            rw_s[1][flag_s["H105r001"]] .= -1.0
            rw_s[2][flag_s["H105r002"]] .= -1.0
        elseif id == "H105r005"
            rw_s[1][flag_s["H105r005"]] .= -1.0
        elseif id == "J303" || id == "J303r003"
            rw_s[1][flag_s["J303r003"]] .= -1.0
        end
        rw = [rw1[k] .* rw2[k] .* rw_s[k] for k =1:length(W)]
        data_r = [data[k] .* rw[k][vcfg[k]] for k=1:length(data)]
        return (data_r, rw)
    end
end

function corr_obs(cdata::CorrData, corr::Corr;
                  real::Bool=true,
                  rw::Union{Array{Float64, 2}, Nothing}=nothing,
                  L::Int64=1, info::Bool=false,
                  idm::Union{Vector{Int64},Nothing}=nothing,
                  nms::Int64=Int64(maximum(cdata.vcfg)),
                  flag_strange::Bool=false)
    real ? data = cdata.re_data ./ L^3 : data = cdata.im_data ./ L^3
    nt = size(data)[2]
    idm = isnothing(idm) ? Int64.(cdata.vcfg) : idm
    if isnothing(rw)
        # idm = isnothing(idm) ? collect(1:nms) : idm
        obs = [uwreal(data[:, x0], cdata.id, idm, nms) for x0 = 1:nt]
    else
        # idm = isnothing(idm) ? collect(1:nms) : idm
        data_r, W = apply_rw(data, rw, cdata.vcfg, id=cdata.id, fs=flag_strange)
        ow = [uwreal(data_r[:, x0], cdata.id, idm, nms) for x0 = 1:nt]
        W_obs = uwreal(W, cdata.id, idm, nms)
        obs = [ow[x0] / W_obs for x0 = 1:nt]
    end
    corr = update(corr,obs=obs)
    if info
        return !isnothing(rw) ?  (corr,obs) : (corr,ow,W_obs)
    else
        return corr
    end
end


    #=
function read_mesons(path::Vector{String}, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false,
    nnoise_trunc::Union{Int64, Nothing}=nothing)
    res = read_mesons.(path, g1, g2, id=id, legacy=legacy, nnoise_trunc=nnoise_trunc)
    nrep = length(res)
    ncorr = length(res[1])

    cdata = Vector{Vector{CData}}(undef, ncorr)
    for icorr = 1:ncorr
        cdata[icorr] = Vector{CData}(undef, nrep)
        for r = 1:nrep
            cdata[icorr][r] = res[r][icorr]
        end
    end
    return cdata
end

@doc raw"""
    read_mesons(path::String, gamma::Vector{Tuple{String,String}}; id::Union{String, Nothing}=nothing, legacy::Bool=false,  nnoise_trunc::Union{Int64, Nothing}=nothing)

    read_mesons(path::Vector{String}, gamma::Vector{Tuple{String,String}}; id::Union{String, Nothing}=nothing, legacy::Bool=false,  nnoise_trunc::Union{Int64, Nothing}=nothing)


Read a meson.dat file and return a `Vector{Cdata}` with for different masses and the  Dirac structures specified in `gamma`.
The Dirac structures are specified into `gamma` in a vector of Tuple as `(Gsrc,Gsnk)`.

This method is best used when muliple, but specific, Dirac structure are requested. It is up to the user reorganized the output in the preferred way.
"""
function read_mesons(path::String, gamma::Vector{Tuple{String,String}}; id::Union{String, Nothing}=nothing, legacy::Bool=false,
  nnoise_trunc::Union{Int64, Nothing}=nothing)

  type = Vector{Tuple{Int64,Int64}}(undef,length(gamma))

  for i in eachindex(gamma)
    G1,G2 = gamma[i]
    type[i] = (findfirst(x-> x==G1, juobs.gamma_name) - 1, findfirst(x-> x==G2, juobs.gamma_name) - 1)
  end


  if isnothing(id)
      bname = basename(path)
      m = findfirst(r"[A-Z][0-9]{3}r[0-9]{3}", bname)
      id = bname[m[1:4]]
  end

  data = open(path, "r")
  g_header = read_GHeader(path)
  c_header = read_CHeader(path, legacy=legacy)

  ncorr = g_header.ncorr
  tvals = g_header.tvals
  nnoise = g_header.nnoise

  nnoise_trunc = isnothing(nnoise_trunc) ? nnoise : min(nnoise, nnoise_trunc)

  fsize = filesize(path)

  datsize = 4 + sum(getfield.(c_header, :dsize)) * tvals * nnoise #data_size / ncnfg
  ncfg = div(fsize - g_header.hsize - sum(getfield.(c_header, :hsize)), datsize) #(total size - header_size) / data_size

  corr_match = findall(x-> (x.type1,x.type2) in type, c_header)

  seek(data, g_header.hsize + sum(c.hsize for c in c_header))

  res = Array{juobs.CData}(undef, length(corr_match))

  data_re = Array{Float64}(undef, length(corr_match), ncfg, tvals)
  data_im = zeros(length(corr_match), ncfg, tvals)
  vcfg = Array{Int32}(undef, ncfg)

  for icfg = 1:ncfg
      vcfg[icfg] = read(data, Int32)
      c=1
      for k = 1:ncorr
          if k in corr_match
              if c_header[k].is_real == 1
                  tmp = Array{Float64}(undef, tvals*nnoise)
                  read!(data, tmp)
                  tmp2 = reshape(tmp, (nnoise, tvals))
                  tmp2 = mean(tmp2[1:nnoise_trunc, :], dims=1)
                  data_re[c, icfg, :] = tmp2[1, :]
              elseif c_header[k].is_real == 0
                  tmp = Array{Float64}(undef, 2*tvals*nnoise)
                  read!(data, tmp)
                  tmp2 = reshape(tmp, (2, nnoise, tvals))
                  tmp2 = mean(tmp2[:, 1:nnoise_trunc, :], dims=2)
                  data_re[c, icfg, :] = tmp2[1, 1, :]
                  data_im[c, icfg, :] = tmp2[2, 1, :]

              end
              c += 1
          else
              seek(data, position(data)  + c_header[k].dsize*tvals*nnoise)
          end


      end
  end
  ## sort by gamma_structure
  idx = Vector{Int64}(undef, length(corr_match))
  is = 1
  for (G1,G2) in type
    aux = findall(x-> (x.type1 == G1 && x.type2==G2),c_header[corr_match])
    if isnothing(aux)
      continue;
    end
    l =length(aux);
    idx[is:is+l-1] = aux
    is += l;
  end

  for i in eachindex(idx)
    c = idx[i]
    res[i] = juobs.CData(c_header[corr_match[c]], vcfg, data_re[c, :, :], data_im[c, :, :], id)
  end
  close(data)

  return res
end

function read_mesons(path::Vector{String}, gamma::Vector{Tuple{String,String}}; id::Union{String, Nothing}=nothing, legacy::Bool=false,
  nnoise_trunc::Union{Int64, Nothing}=nothing)
  res = [read_mesons(p, gamma, id=id, legacy=legacy, nnoise_trunc=nnoise_trunc) for p in path]
  nrep = length(res)
  ncorr = length(res[1])

  cdata = Vector{Vector{juobs.CData}}(undef, ncorr)
  for icorr = 1:ncorr
      cdata[icorr] = Vector{juobs.CData}(undef, nrep)
      for r = 1:nrep
          cdata[icorr][r] = res[r][icorr]
      end
  end
  return cdata
end



function read_mesons_correction(path::String, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false,
    nnoise_trunc::Union{Int64, Nothing}=nothing)
    t1 = isnothing(g1) ? nothing : findfirst(x-> x==g1, gamma_name) - 1
    t2 = isnothing(g2) ? nothing : findfirst(x-> x==g2, gamma_name) - 1
    if isnothing(id)
        bname = basename(path)
        m = findfirst(r"[A-Z][0-9]{3}r[0-9]{3}", bname)
        id = bname[m[1:4]]
        #id = parse(Int64, bname[m[2:4]])
    end

    data = open(path, "r")
    g_header = read_GHeader(path)
    c_header = read_CHeader(path, legacy=legacy)

    ncorr = g_header.ncorr
    tvals = g_header.tvals
    nnoise = g_header.nnoise

    nnoise_trunc = isnothing(nnoise_trunc) ? nnoise : min(nnoise, nnoise_trunc)

    fsize = filesize(path)

    datsize = 4 + sum(getfield.(c_header, :dsize)) * tvals * nnoise #data_size / ncnfg
    ncfg = div(fsize - g_header.hsize - sum(getfield.(c_header, :hsize)), datsize) #(total size - header_size) / data_size

    corr_match = findall(x-> (x.type1==t1 || isnothing(t1)) && (x.type2==t2 || isnothing(t2)), c_header)


    seek(data, g_header.hsize + sum(c.hsize for c in c_header))

    res = Array{CData}(undef, div(length(corr_match), 2)) # Modification: total length is divided by 2

    data_re = zeros(div(length(corr_match), 2), ncfg, tvals) # Modification: total length is divided by 2
    data_im = zeros(div(length(corr_match), 2), ncfg, tvals) # Modification: total length is divided by 2
    vcfg = Array{Int32}(undef, ncfg)

    for icfg = 1:ncfg
        vcfg[icfg] = read(data, Int32)
        c = 1
        sgn = +1 # sign it changes at ncorr / 2. O_exact - O_sloppy
        for k = 1:ncorr
            if k in corr_match
                if c_header[k].is_real == 1
                    tmp = Array{Float64}(undef, tvals*nnoise)
                    read!(data, tmp)
                    tmp2 = reshape(tmp, (nnoise, tvals))
                    tmp2 = mean(tmp2[1:nnoise_trunc, :], dims=1)

                    data_re[c, icfg, :] = data_re[c, icfg, :] + sgn * tmp2[1, :]
                elseif c_header[k].is_real == 0
                    tmp = Array{Float64}(undef, 2*tvals*nnoise)
                    read!(data, tmp)
                    tmp2 = reshape(tmp, (2, nnoise, tvals))
                    tmp2 = mean(tmp2[:, 1:nnoise_trunc, :], dims=2)
                    data_re[c, icfg, :] = data_re[c, icfg, :] + sgn * tmp2[1, 1, :]
                    data_im[c, icfg, :] = data_im[c, icfg, :] + sgn * tmp2[2, 1, :]

                end
                c += 1
            else
                seek(data, position(data)  + c_header[k].dsize*tvals*nnoise)
            end
            if k == div(ncorr, 2)
                c = 1
                sgn = -1
            end
        end
    end


    for c in 1:div(length(corr_match), 2)
    res[c] = juobs.CData(c_header[corr_match[c]], vcfg, data_re[c, :, :], data_im[c, :, :], id)
  end

    close(data)

    return res
end


@doc raw"""
    read_mesons_correction(path::String, gamma::Vector{Tuple{String,String}}; id::Union{String, Nothing}=nothing, legacy::Bool=false,  nnoise_trunc::Union{Int64, Nothing}=nothing)

    read_mesons_correction(path::Vector{String}, gamma::Vector{Tuple{String,String}}; id::Union{String, Nothing}=nothing, legacy::Bool=false,  nnoise_trunc::Union{Int64, Nothing}=nothing)


Read a meson.dat file with correction data and return a `Vector{Cdata}` with for different masses and the  Dirac structures specified in `gamma`. (See also [`read_meson`](@ref))
The Dirac structures are specified into `gamma` in a vector of Tuple as `(Gsrc,Gsnk)`.

This method is best used when muliple, but specific, Dirac structure are requested. It is up to the user reorganized the output in the preferred way.
"""
function read_mesons_correction(path::Vector{String}, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false,
    nnoise_trunc::Union{Int64, Nothing}=nothing)
    res = read_mesons_correction.(path, g1, g2, id=id, legacy=legacy, nnoise_trunc=nnoise_trunc)
    nrep = length(res)
    ncorr = length(res[1])

    cdata = Vector{Vector{CData}}(undef, ncorr)
    for icorr = 1:ncorr
        cdata[icorr] = Vector{CData}(undef, nrep)
        for r = 1:nrep
            cdata[icorr][r] = res[r][icorr]
        end
    end
    return cdata
end

function read_mesons_correction(path::String, gamma::Vector{Tuple{String,String}}; id::Union{String, Nothing}=nothing, legacy::Bool=false,
  nnoise_trunc::Union{Int64, Nothing}=nothing)

  type = Vector{Tuple{Int64,Int64}}(undef,length(gamma))

  for i in eachindex(gamma)
    G1,G2 = gamma[i]
    type[i] = (findfirst(x-> x==G1, juobs.gamma_name) - 1, findfirst(x-> x==G2, juobs.gamma_name) - 1)
  end

  if isnothing(id)
      bname = basename(path)
      m = findfirst(r"[A-Z][0-9]{3}r[0-9]{3}", bname)
      id = bname[m[1:4]]
      #id = parse(Int64, bname[m[2:4]])
  end

  data = open(path, "r")
  g_header = juobs.read_GHeader(path)
  c_header = juobs.read_CHeader(path, legacy=legacy)

  ncorr = g_header.ncorr
  tvals = g_header.tvals
  nnoise = g_header.nnoise

  nnoise_trunc = isnothing(nnoise_trunc) ? nnoise : min(nnoise, nnoise_trunc)

  fsize = filesize(path)

  datsize = 4 + sum(getfield.(c_header, :dsize)) * tvals * nnoise #data_size / ncnfg
  ncfg = div(fsize - g_header.hsize - sum(getfield.(c_header, :hsize)), datsize) #(total size - header_size) / data_size

  corr_match = findall(x-> (x.type1,x.type2) in type, c_header)


  seek(data, g_header.hsize + sum(c.hsize for c in c_header))

  res = Array{juobs.CData}(undef, div(length(corr_match), 2)) # Modification: total length is divided by 2

  data_re = zeros(div(length(corr_match), 2), ncfg, tvals) # Modification: total length is divided by 2
  data_im = zeros(div(length(corr_match), 2), ncfg, tvals) # Modification: total length is divided by 2
  vcfg = Array{Int32}(undef, ncfg)

  for icfg = 1:ncfg
      vcfg[icfg] = read(data, Int32)
      c = 1
      sgn = +1 # sign it changes at ncorr / 2. O_exact - O_sloppy
      for k = 1:ncorr
          if k in corr_match
              if c_header[k].is_real == 1
                  tmp = Array{Float64}(undef, tvals*nnoise)
                  read!(data, tmp)
                  tmp2 = reshape(tmp, (nnoise, tvals))
                  tmp2 = mean(tmp2[1:nnoise_trunc, :], dims=1)

                  data_re[c, icfg, :] = data_re[c, icfg, :] + sgn * tmp2[1, :]
              elseif c_header[k].is_real == 0
                  tmp = Array{Float64}(undef, 2*tvals*nnoise)
                  read!(data, tmp)
                  tmp2 = reshape(tmp, (2, nnoise, tvals))
                  tmp2 = mean(tmp2[:, 1:nnoise_trunc, :], dims=2)
                  data_re[c, icfg, :] = data_re[c, icfg, :] + sgn * tmp2[1, 1, :]
                  data_im[c, icfg, :] = data_im[c, icfg, :] + sgn * tmp2[2, 1, :]

              end
              c += 1
          else
              seek(data, position(data)  + c_header[k].dsize*tvals*nnoise)
          end
          if k == div(ncorr, 2)
              c = 1
              sgn = -1
          end
      end
  end
  for c = 1:div(length(corr_match), 2)
      res[c] = juobs.CData(c_header[corr_match[c]], vcfg, data_re[c, :, :], data_im[c, :, :], id)
  end
  close(data)

  return res
end

function read_mesons_correction(path::Vector{String}, gamma::Vector{Tuple{String,String}}; id::Union{String, Nothing}=nothing, legacy::Bool=false,
nnoise_trunc::Union{Int64, Nothing}=nothing)
res = [read_mesons_correction(p, gamma, id=id, legacy=legacy, nnoise_trunc=nnoise_trunc) for p in path]
nrep = length(res)
ncorr = length(res[1])

cdata = Vector{Vector{juobs.CData}}(undef, ncorr)
for icorr = 1:ncorr
    cdata[icorr] = Vector{juobs.CData}(undef, nrep)
    for r = 1:nrep
        cdata[icorr][r] = res[r][icorr]
    end
end
return cdata
end

@doc raw"""
    read_mesons_chunks(paths::Vector{String}, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false,
    nnoise_trunc::Union{Int64, Nothing}=nothing)

This function reads mesons dat files stored in chunks and returns a single vector of `CData` structures for different masses and Dirac structures.
  the paths to each chunks is given through the variable `paths::Vector{String}`
  Dirac structures `g1` and/or `g2` can be passed as string arguments in order to filter correlators.
  ADerrors id can be specified as argument. If is not specified, the `id` is fixed according to the ensemble name (example: "H400"-> id = "H400")

  *For the old version (without smearing, distance preconditioning and theta) set legacy=true.

  Examples:
  ```@example
  read_mesons([path_chunk1,path_chunk2,path_chunk3])
  read_mesons([path_chunk1,path_chunk2,path_chunk3], "G5")
  read_mesons([path_chunk1,path_chunk2,path_chunk3], nothing, "G5")
  read_mesons([path_chunk1,path_chunk2,path_chunk3], "G5", "G5")
  read_mesons([path_chunk1,path_chunk2,path_chunk3], "G5", "G5", id="H100")
  read_mesons([path_chunk1,path_chunk2,path_chunk3], "G5_d2", "G5_d2", legacy=true)
  ```
"""
function read_mesons_chunks(paths::Vector{String}, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false,
  nnoise_trunc::Union{Int64, Nothing}=nothing)
  data = read_mesons(paths[1],g1,g2,id=id,legacy=legacy, nnoise_trunc=nnoise_trunc)
  for i in 2:length(paths)
    temp = read_mesons(paths[i],g1,g2 ,id=id,legacy=legacy, nnoise_trunc=nnoise_trunc)
    concat_data!(data,temp)
  end
  return data
end

function read_rw(path::String; v::String="1.2")
    data = open(path, "r")
    nrw = read(data, Int32)

    if v == "1.4" || v =="1.6"
        nfct = Array{Int32}(undef, nrw)
        read!(data, nfct)
        nfct_inheader = 1
    elseif v=="1.2"
        nfct = ones(Int32, nrw)
        nfct_inheader = 0
    else
        error("Version not supported")
    end
    nsrc = Array{Int32}(undef, nrw)
    read!(data, nsrc)
    glob_head_size = 4 + 4*nrw*(nfct_inheader + 1)
    datsize = 4 + 2*8*sum(nsrc .* nfct)

    fsize = filesize(path)

    ncnfg = Int32((fsize - glob_head_size)/datsize)
    r_data = Array{Array{Float64}}(undef, nrw)

    [r_data[k] = zeros(Float64, nfct[k], nsrc[k], ncnfg) for k = 1:nrw]
    vcfg = Array{Int32}(undef, ncnfg)
    for icfg = 1:ncnfg
        vcfg[icfg] = read(data, Int32)
        for irw = 1:nrw
            for ifct = 1:nfct[irw]
                tmp = zeros(Float64, nsrc[irw])
                seek(data, position(data) + 8 * nsrc[irw])
                read!(data, tmp)
                r_data[irw][ifct, :, icfg] = tmp[:]
            end
        end
    end
    close(data)
    return r_data
end

@doc raw"""
read_rw_openQCD2(path::String; print_info::Bool=false)

This function reads the reweighting factors generated with openQCD version 2.#.
The flag print_info if set to true print additional information for debugging
"""
function read_rw_openQCD2(path::String; print_info::Bool=false)

    data = open(path, "r")
    nrw = read(data, Int32)
    nrw = Int32(nrw / 2)

    nfct = Array{Int32}(undef, nrw)
    read!(data, nfct)

    nsrc = Array{Int32}(undef, nrw)
    read!(data, nsrc)
    null = read(data, Int32)
    if null !== Int32(0)
        error("In openQCD 2.0 this Int32 should be a zero.")
    end

    data_array = Array{Array{Float64}}(undef, nrw)
    [data_array[k] = Array{Float64}(undef, 0) for k in 1:nrw]
    vcfg = Vector{Int32}(undef, 0)
    while !eof(data)

        push!(vcfg, read(data, Int32))
        if print_info
            println("\n ######## cnfg: ", vcfg[end])
        end

        for k in 1:nrw
            read_array_rwf_dat_openQCD2(data)
            tmp_rw, n = read_array_rwf_dat_openQCD2(data)

            tmp_nfct=1.0
            for j in 1:n[1]
                tmp_nfct *= mean((exp.(.-tmp_rw[j])))
            end
            push!(data_array[k], tmp_nfct)
        end
    end

    return permutedims(hcat(data_array...), (2,1))
end


function read_array_rwf_dat_openQCD2(data::IOStream; print_info::Bool=false)

    d = read(data, Int32)
    n = Array{Int32}(undef, d)
    read!(data, n)
    size = read(data, Int32)
    m = prod(n)

    if print_info
        println("d: ", d)
        println("n: ", n)
        println("size: ", size)
        println("m: ", m)
    end

    if size == 4
        types = Int32
    elseif size == 8
        types = Float64
    elseif size == 16
        types = Float64x2
    else
        error("No type with size=$(size) supported")
    end

    tmp_data = Array{types}(undef, m)
    read!(data, tmp_data)

    res = parse_array_openQCD2(d, n, tmp_data, quadprec=true)

    return res, n
end

function parse_array_openQCD2(d, n, dat; quadprec=true)

    if d != 2
        error("dat must be a two-dimensional array")
    end
    res = Vector{Vector{Float64}}(undef, 0)

    for k in range(1,n[1])
        tmp = dat[(k-1)*n[2]+1:k*n[2]]
        if quadprec
            tmp2 = Vector{Float64}(undef, 0)
            for j in range(start=1,step=2,stop=length(tmp))
                push!(tmp2, tmp[j])
            end
            push!(res, tmp2)
        else
            push!(res, tmp)
        end
    end

    return res
end

@doc raw"""
    read_ms1(path::String; v::String="1.2")

Reads openQCD ms1 dat files at a given path. This method returns a matrix `W[irw, icfg]` that contains the reweighting factors, where
`irw` is the `rwf` index and icfg the configuration number.
The function is compatible with the output files of openQCD v=1.2, 1.4, 1.6 and 2.0. Version can be specified as argument.

Examples:
```@example
read_ms1(path)
read_ms1(path, v="1.4")
read_ms1(path, v="1.6")
read_ms1(path, v="2.0")
```
"""
function read_ms1(path::String; v::String="1.2")

    if v == "2.0"
        return read_rw_openQCD2(path)
    end
    r_data = read_rw(path, v=v)
    nrw = length(r_data)
    ncnfg = size(r_data[1])[3]
    W = zeros(Float64, nrw, ncnfg)
    [W[k, :] = prod(mean(exp.(.-r_data[k]), dims=2), dims=1) for k = 1:nrw]
    return W
end
@doc raw"""
    read_md(path::String)

Reads openQCD  pbp.dat files at a given path. This method returns a matrix `md[irw, icfg]` that contains the derivatives ``dS/dm``, where
``md[irw=1] = dS/dm_l`` and ``md[irw=2] = dS/dm_s``

``Seff = -tr(log(D+m))``

``dSeff/ dm = -tr((D+m)^-1)``

Examples:
```@example
md = read_md(path)
```
"""
function read_md(path::String)
    r_data = read_rw(path, v="1.4")
    nrw = length(r_data)
    ncnfg = size(r_data[1])[3]
    md = zeros(Float64, nrw, ncnfg)
    [md[k, :] = prod(.-mean(r_data[k], dims=2), dims=1) for k = 1:nrw]
    return md
end

@doc raw"""
    read_ms(path::String; id::Union{String, Nothing}=nothing, dtr::Int64=1, obs::String="Y")

Reads openQCD ms dat files at a given path. This method return YData:

- `t(t)`: flow time values

- `obs(icfg, x0, t)`: the time-slice sums of the densities of the observable (Wsl, Ysl or Qsl)

- `vtr`: vector that contains trajectory number

- `id`: ensmble id

`dtr` = `dtr_cnfg` / `dtr_ms`, where `dtr_cnfg` is the number of trajectories computed before saving the configuration. `dtr_ms`
is the same but applied to the ms.dat file.

Examples:
```@example
Y = read_ms(path)
```
"""
function read_ms(path::String; id::Union{String, Nothing}=nothing, dtr::Int64=1 , obs::String="Y")
    if isnothing(id)
        bname = basename(path)
        m = findfirst(r"[A-Z][0-9]{3}r[0-9]{3}", bname)
        id = bname[m[1:4]]
        #id = parse(Int64, bname[m[2:4]])
    end
    data = open(path, "r")

    dn = read(data, Int32)
    nn = read(data, Int32)
    tvals = read(data, Int32)
    eps = read(data, Float64)

    fsize = filesize(path)
    datsize=4 + 3*8*(nn + 1) * tvals # measurement size of each cnfg

    ntr = Int32((fsize - 3*4 - 8) / datsize)

    if mod(ntr, dtr) != 0
        println("ntr = ", ntr)
        println("dtr = ", dtr)
        @warn("ntr / dtr must be exact, truncating...")
    end

    vntr = Vector{Int32}(undef, div(ntr, dtr))

    # x0, t, cfg
    Wsl = Array{Float64}(undef, div(ntr, dtr), tvals, nn + 1)
    Ysl = Array{Float64}(undef, div(ntr, dtr), tvals, nn + 1)
    Qsl = Array{Float64}(undef, div(ntr, dtr), tvals, nn + 1)

    k = 0
    for itr = 1:ntr
        tmp = read(data, Int32)
        if mod(itr, dtr) == 0
            k += 1
            vntr[k] = tmp
        end

        for iobs = 1:3
            for inn = 0:nn
                tmp2 = Vector{Float64}(undef, tvals)
                read!(data, tmp2)
                if mod(itr, dtr) == 0
                    if iobs == 1
                        Wsl[k, :, inn + 1] = tmp2
                    elseif iobs == 2
                        Ysl[k, :, inn + 1] = tmp2
                    elseif iobs == 3
                        Qsl[k, :, inn + 1] = tmp2
                    end
                end
            end
        end
    end
    close(data)
    t = Float64.(0:nn) .* dn .* eps

    if obs == "W"
        return YData(vntr, t, Wsl, id)
    elseif obs == "Y"
        return YData(vntr, t, Ysl, id)
    elseif obs == "Q"
        return YData(vntr, t, Qsl, id)
    else
        println("obs = ", obs," is not valid")
        return nothing
    end
end

@doc raw"""
    truncate_data!(data::YData, nc::Int64)

    truncate_data!(data::Vector{YData}, nc::Vector{Int64})

    truncate_data!(data::Vector{CData}, nc::Int64)

    truncate_data!(data::Vector{Vector{CData}}, nc::Vector{Int64})

Truncates the output of `read_mesons` and `read_ms` taking the first `nc` configurations.

Examples:
```@example
#Single replica
dat = read_mesons(path, "G5", "G5")
Y = read_ms(path)
truncate_data!(dat, nc)
truncate_data!(Y, nc)

#Two replicas
dat = read_mesons([path1, path2], "G5", "G5")
Y = read_ms.([path1_ms, path2_ms])
truncate_data!(dat, [nc1, nc2])
truncate_data!(Y, [nc1, nc2])
```
"""
function truncate_data!(data::YData, nc::Int64)
    data.vtr = data.vtr[1:nc]
    data.obs = data.obs[1:nc, :, :]
    return nothing
end
function truncate_data!(data::Vector{YData}, nc::Vector{Int64})
    truncate_data!.(data, nc)
    return nothing
end
function truncate_data!(data::Vector{CData}, nc::Int64)
    N = length(data)
    for k = 1:N
        data[k].vcfg = data[k].vcfg[1:nc]
        data[k].re_data = data[k].re_data[1:nc, :]
        data[k].im_data = data[k].im_data[1:nc, :]
    end
    return nothing
end
function truncate_data!(data::Vector{Vector{CData}}, nc::Vector{Int64})
    N = length(data)
    R = length(data[1])
    for k = 1:N
        for r = 1:R
            data[k][r].vcfg = data[k][r].vcfg[1:nc[r]]
            data[k][r].re_data = data[k][r].re_data[1:nc[r], :]
            data[k][r].im_data = data[k][r].im_data[1:nc[r], :]
        end
    end
    return nothing
end

@doc raw"""

    concat_data!(data1::Vector{CData}, data2::Vector{CData})

    concat_data!(data1::Vector{Vector{CData}}, data2::Vector{Vector{CData}})

Concatenates the output of 2 calls to `read_mesons` over configurations.
Both data must have the same number of replicas and correlators.
The output is saved in the first argument, so if you want to concatenate
3 data sets: `concat_data!(data1, data2); concat_data!(data1, data3)`

Examples:
```@example
#Single replica
dat = read_mesons(path, "G5", "G5")
dat2 = read_mesons(path2, "G5", "G5")
concat_data!(dat, dat2)

#Two replicas
dat = read_mesons([path1, path2], "G5", "G5")
dat2 = read_mesons([path3, path4], "G5", "G5")
concat_data!(dat, dat2)
```
"""
function concat_data!(data1::Vector{juobs.CData}, data2::Vector{juobs.CData})
    N = length(data1)
    if length(data1) != length(data2)
        error("number of correlators do not match")
    end
    for k = 1:N
        data1[k].vcfg = vcat(data1[k].vcfg, data2[k].vcfg)
        data1[k].re_data = vcat(data1[k].re_data, data2[k].re_data)
        data1[k].im_data = vcat(data1[k].im_data, data2[k].im_data)
        idx = sortperm(data1[k].vcfg)
        data1[k].vcfg = data1[k].vcfg[idx]
        data1[k].re_data = data1[k].re_data[idx, :]
        data1[k].im_data = data1[k].im_data[idx, :]
    end
    return nothing
end

function concat_data!(data1::Vector{Vector{juobs.CData}}, data2::Vector{Vector{juobs.CData}})
    N = length(data1)
    if length(data1) != length(data2)
        error("number of correlators do not match")
    end
    R = length(data1[1])
    if length(data1[1]) != length(data2[1])
        error("number of replicas do not match")
    end
    for k = 1:N
        for r = 1:R
            data1[k][r].vcfg = vcat(data1[k][r].vcfg, data2[k][r].vcfg)
            data1[k][r].re_data = vcat(data1[k][r].re_data, data2[k][r].re_data)
            data1[k][r].im_data = vcat(data1[k][r].im_data, data2[k][r].im_data)
            idx = sortperm(data1[k][r].vcfg)
            data1[k][r].vcfg = data1[k][r].vcfg[idx]
            data1[k][r].re_data = data1[k][r].re_data[idx, :]
            data1[k][r].im_data = data1[k][r].im_data[idx, :]

        end
    end
    return nothing
end

function get_YM(path::String, ens::EnsInfo; rw=false, ws::ADerrors.wspace=ADerrors.wsg, w0::Union{Float64, Nothing}=nothing)

    path_ms = joinpath(path, ens.id, "gf")
    path_ms = filter(x->occursin(".dat", x), readdir(path_ms, join=true))
    Y = read_ms.(path_ms, dtr=ens.dtr)

    nr = length(Y)
    Ysl = getfield.(Y, :obs)
    t = getfield.(Y, :t)
    t = t[1]
    id = getfield.(Y, :id)
    replica = size.(Ysl, 1)

    L = ens.L
    id = ens.id
    #T = length(Y[:,1]) - y0
    y0 = 1 ## assumes this is the case, hardcoded, some ensembles will not fulfil !
    println("WARNING!: make sure t_src is 1 in this ensemble")

    #Truncation
    if id in keys(ADerrors.wsg.str2id)
        n_ws = findfirst(x-> x == ws.str2id[id], ws.map_nob)
        if !isnothing(n_ws)
            ivrep_ws = ws.fluc[n_ws].ivrep

            if length(replica) != length(ivrep_ws)
                error("Different number of replicas")
            end

            for k = 1:length(replica)
                if replica[k] > ivrep_ws[k]
                    println("Automatic truncation in Ysl ", ivrep_ws[k], " / ", replica[k], ". R = ", k)
                    Ysl[k] = Ysl[k][1:ivrep_ws[k], :, :]
                elseif replica[k] < ivrep_ws[k]
                    error("Automatic truncation failed. R = ", replica[k], "\nTry using truncate_data!")
                end
            end
            replica = size.(Ysl, 1)
        end
    end

    tmp = Ysl[1]
    [tmp = cat(tmp, Ysl[k], dims=1) for k = 2:nr]
    xmax = size(tmp, 2)
    T = xmax - 1 - y0

    Y_aux = Matrix{uwreal}(undef, xmax, length(t))

    if rw
        path_rw = joinpath(path, ens.id, "rwf")
        path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
        if ens.id == "D200"
            rwf_1 = read_ms1.([path_rw[1]], v=ens.vrw)
            rwf_2 = read_ms1.([path_rw[2]], v=ens.vrw)
            rwf = [hcat(rwf_2[1],rwf_1[1])]
        elseif ens.id in ["E250", "E300", "J500", "J501", "D450"]
            rwf = [read_ms1(path_rw[i], v=ens.vrw[i]) for i in 1:length(ens.vrw)]
            [Ysl[k] = Ysl[k][1:size(rwf[k],2), :, :] for k in 1:length(Ysl)]
        else
            rwf = read_ms1.(path_rw, v=ens.vrw)
        end

        Ysl_r, W = juobs.apply_rw(Ysl, rwf)
        tmp_r = Ysl_r[1]
        tmp_W = W[1]
        [tmp_r = cat(tmp_r, Ysl_r[k], dims=1) for k = 2:nr]
        [tmp_W = cat(tmp_W, W[k], dims=1) for k = 2:nr]
        W_obs = uwreal(tmp_W, id, replica, collect(1:length(tmp_W)), sum(replica))
        WY_aux = Matrix{uwreal}(undef, xmax, length(t))
    end
    for i = 1:xmax
        k = 1
        for j = 1:length(t)
            if !rw
                Y_aux[i, k] = uwreal(tmp[:, i, j], id, replica)
            else
                WY_aux[i, k] = uwreal(tmp_r[:, i, j], id, replica, collect(1:length(tmp_W)), sum(replica))
                Y_aux[i, k] = WY_aux[i, k] / W_obs
            end
            k = k + 1
        end
    end

    t2YM = similar(Y_aux)
    tdt2YM = similar(Y_aux)
    for i in 1:length(Y_aux[:,1])
        t2YM[i,:] = Y_aux[i,:] .* t .^ 2 ./ L ^ 3
    end
    for i in 1:length(t2YM[:,1])
        if isnothing(w0)
            tdt2YM[i,2:end-1] = [(t2YM[i,j+1] - t2YM[i,j-1]) / (t[j+1] - t[j-1]) * t[j] for j in 2:length(t2YM[i,:])-1]
            tdt2YM[i,1] = tdt2YM[i,end] = t2YM[i,1]
        else
            ixm = findmin(abs.(t .- (w0-0.5)))[2]
            ixM = findmin(abs.(t[1:end-1] .- (w0+0.5)))[2]
            tdt2YM[i,1:end] .= t2YM[i,2]
            tdt2YM[i,ixm:ixM] = [(t2YM[i,j+1] - t2YM[i,j-1]) / (t[j+1] - t[j-1]) * t[j] for j in ixm:ixM]
        end
    end

    return t2YM, tdt2YM, W_obs, t
end

function get_YM_dYM(path::String, ens::EnsInfo; rw=false, ws::ADerrors.wspace=ADerrors.wsg, w0::Union{Float64, Nothing}=nothing, tau::Union{Float64, Nothing}=nothing)

    path_ms = joinpath(path, ens.id, "gf")
    path_ms = filter(x->occursin(".dat", x), readdir(path_ms, join=true))
    Y = read_ms.(path_ms, dtr=ens.dtr)

    truncate_data!(Y, ens.cnfg)

    nr = length(Y)
    Ysl = getfield.(Y, :obs)
    t = getfield.(Y, :t)
    t = t[1]
    id = getfield.(Y, :id)
    replica = size.(Ysl, 1)

    L = ens.L
    id = ens.id
    #T = length(Y[:,1]) - y0
    y0 = 1 ## assumes this is the case, hardcoded, some ensembles will not fulfil !
    println("WARNING!: make sure t_src is 1 in this ensemble")

    #Truncation
    if id in keys(ADerrors.wsg.str2id)
        n_ws = findfirst(x-> x == ws.str2id[id], ws.map_nob)
        if !isnothing(n_ws)
            ivrep_ws = ws.fluc[n_ws].ivrep

            if length(replica) != length(ivrep_ws)
                error("Different number of replicas")
            end

            for k = 1:length(replica)
                if replica[k] > ivrep_ws[k]
                    println("Automatic truncation in Ysl ", ivrep_ws[k], " / ", replica[k], ". R = ", k)
                    Ysl[k] = Ysl[k][1:ivrep_ws[k], :, :]
                elseif replica[k] < ivrep_ws[k]
                    error("Automatic truncation failed. R = ", replica[k], "\nTry using truncate_data!")
                end
            end
            replica = size.(Ysl, 1)
        end
    end

    tmp = Ysl[1]
    [tmp = cat(tmp, Ysl[k], dims=1) for k = 2:nr]
    xmax = size(tmp, 2)
    T = xmax - 1 - y0

    Y_aux = Matrix{uwreal}(undef, xmax, length(t))

    if rw
        path_rw = joinpath(path, ens.id, "rwf")
        path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
        if ens.id in ["H105", "J500", "J501"]
            rwf = [read_ms1(path_rw[i], v=ens.vrw[i]) for i in 1:length(ens.vrw)]
            [Ysl[k] = Ysl[k][1:size(rwf[k],2), :, :] for k in 1:length(Ysl)]
        else
            rwf = read_ms1.(path_rw, v=ens.vrw)
        end

        Ysl_r, W = juobs.apply_rw(Ysl, rwf)
        tmp_r = Ysl_r[1]
        tmp_W = W[1]
        [tmp_r = cat(tmp_r, Ysl_r[k], dims=1) for k = 2:nr]
        [tmp_W = cat(tmp_W, W[k], dims=1) for k = 2:nr]
        W_obs = uwreal(tmp_W, id, replica, collect(1:length(tmp_W)), sum(replica))
        WY_aux = Matrix{uwreal}(undef, xmax, length(t))
    end
    for i = 1:xmax
        k = 1
        for j = 1:length(t)
            if !rw
                Y_aux[i, k] = uwreal(tmp[:, i, j], id, replica)
            else
                WY_aux[i, k] = uwreal(tmp_r[:, i, j], id, replica, collect(1:length(tmp_W)), sum(replica))
                Y_aux[i, k] = WY_aux[i, k] / W_obs
            end
            k = k + 1
        end
    end

    if tau == nothing
        t2YM = similar(Y_aux)
        tdt2YM = Matrix{uwreal}(undef,size(t2YM)[1],size(t2YM)[2]-2)
    else
        tau_ind = Int64(div(tau, t[2]-t[1]))
        if tau >= 0.0
            global t = t[1:end-tau_ind]
        elseif tau < 0.0
            global t = t[1-tau_ind:end]
        end
        t2YM = Matrix{uwreal}(undef,size(Y_aux)[1],length(t))
        tdt2YM = Matrix{uwreal}(undef,size(t2YM)[1],size(t2YM)[2]-2)
    end
    for i in 1:length(Y_aux[:,1])
        if tau == nothing
            t2YM[i,:] = Y_aux[i,:] .* t .^ 2 ./ L ^ 3
        else
            if tau >= 0.0
                t2YM[i,1:end] = Y_aux[i,1+tau_ind:end] .* t .^ 2 ./ L ^ 3
            elseif tau < 0.0
                t2YM[i,1:end] = Y_aux[i,1:end+tau_ind] .* t .^ 2 ./ L ^ 3
            end
        end
    end
    for i in 1:length(t2YM[:,1])
        if isnothing(w0)
            tdt2YM[i,1:end] = [(t2YM[i,j+1] - t2YM[i,j-1]) / (t[j+1] - t[j-1]) * t[j] for j in 2:length(t2YM[i,:])-1]
            global t_aux = t[2:end-1]
        else
            global t_aux = t[2:end-1]
            ixm = findmin(abs.(t_aux .- (w0-0.5)))[2]
            ixM = findmin(abs.(t_aux[1:end-1] .- (w0+0.5)))[2]
            tdt2YM[i,1:ixm-1] .= t2YM[i,2]
            tdt2YM[i,ixm:end] = [(t2YM[i,j+1] - t2YM[i,j-1]) / (t[j+1] - t[j-1]) * t[j] for j in ixm+1:length(t2YM[i,:])-1]
        end
    end

    return t2YM, tdt2YM, W_obs, t, t_aux
end

function read_mesons_multichunks(path::String, ens::String, g1::String="G5", g2::String="G5"; legacy::Bool=false)
    idx_ts001 = findall(x->occursin("ts001",x), db[ens])
    idx_tsT = findall(x->!occursin("ts001",x), db[ens])

    store_cdata_aux = []
    for i in db[ens]
        aux = filter(x->occursin(i, x), readdir(path, join=true))
        data_chunk = read_mesons(aux, g1, g2, legacy=legacy);
        push!(store_cdata_aux, data_chunk)
    end

    #=
    if ens == "E300"
        for i in 1:length(idx_ts001)-1
            concat_data!(store_cdata_aux[idx_ts001[1]][1:2:end], store_cdata_aux[idx_ts001[i+1]])
        end
        concat_data!(store_cdata_aux[idx_tsT[1]], store_cdata_aux[idx_ts001[1]][2:2:end])
        for i in 1:length(idx_tsT)-1
            concat_data!(store_cdata_aux[idx_tsT[1]], store_cdata_aux[idx_tsT[i+1]])
        end
        store_cdata_aux[idx_ts001[1]] = store_cdata_aux[idx_ts001[1]][1:2:end]
    else
        for i in 1:length(idx_ts001)-1
            concat_data!(store_cdata_aux[idx_ts001[1]], store_cdata_aux[idx_ts001[i+1]])
        end
        for i in 1:length(idx_tsT)-1
            concat_data!(store_cdata_aux[idx_tsT[1]], store_cdata_aux[idx_tsT[i+1]])
        end
    end
    =#
    for i in 1:length(idx_ts001)-1
        concat_data!(store_cdata_aux[idx_ts001[1]], store_cdata_aux[idx_ts001[i+1]])
    end
    for i in 1:length(idx_tsT)-1
        concat_data!(store_cdata_aux[idx_tsT[1]], store_cdata_aux[idx_tsT[i+1]])
    end
    dat_ts001 = store_cdata_aux[idx_ts001[1]]
    dat_ts190 = store_cdata_aux[idx_tsT[1]]

    dat = Vector{Vector{juobs.CData}}()
    for i in 1:length(dat_ts001)
        push!(dat, dat_ts001[i])
        push!(dat, dat_ts190[i])
    end

    return dat
end

function read_mesons_correction_multichunks(path::String, ens::String, g1::String="G5", g2::String="G5"; legacy::Bool=false)
    idx_ts001 = findall(x->occursin("ts001",x), db_c[ens])
    idx_tsT = findall(x->!occursin("ts001",x), db_c[ens])

    store_cdata_aux = []
    for i in db_c[ens]
        aux = filter(x->occursin(i, x), readdir(path, join=true))
        data_chunk = read_mesons_correction(aux, g1, g2, legacy=legacy);
        push!(store_cdata_aux, data_chunk)
    end

    for i in 1:length(idx_ts001)-1
        concat_data!(store_cdata_aux[idx_ts001[1]], store_cdata_aux[idx_ts001[i+1]])
    end
    for i in 1:length(idx_tsT)-1
        concat_data!(store_cdata_aux[idx_tsT[1]], store_cdata_aux[idx_tsT[i+1]])
    end
    dat_ts001 = store_cdata_aux[idx_ts001[1]]
    if length(idx_tsT) >= 1
        dat_ts190 = store_cdata_aux[idx_tsT[1]]
        dat = Vector{Vector{juobs.CData}}()
        for i in 1:length(dat_ts001)
            push!(dat, dat_ts001[i])
            push!(dat, dat_ts190[i])
        end
    else
        dat = dat_ts001
    end

    return dat
end

function get_corr_TSM_multichunks(path::String, ens::EnsInfo; info=false)
    path = joinpath(path, ens.id)
    path_sl = joinpath.(path, "sloppy")
    path_c = joinpath.(path, "correc")
    path_rw = joinpath(path, "rwf_def")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    if length(path_rw) == 0
        path_rw = joinpath(path, "rwf")
        global path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    end

    pp_dat = read_mesons_multichunks(path_sl, ens.id, "G5", "G5")
    ap_dat = read_mesons_multichunks(path_sl, ens.id, "G5", "G0G5")
    pp_dat_c = read_mesons_correction_multichunks(path_c, ens.id, "G5", "G5")
    ap_dat_c = read_mesons_correction_multichunks(path_c, ens.id, "G5", "G0G5")

    rwf = [read_ms1(path_rw[i], v=ens.vrw[i]) for i in 1:length(ens.vrw)]
    cnfg_rw = size.(rwf,2)
    cnfg_trunc_ts001 = [findall(pp_dat[1][i].vcfg .< cnfg_rw[i])[end] for i in 1:length(cnfg_rw)] ## rwf missing some configs at the end
    cnfg_trunc_ts001_c = [findall(pp_dat_c[1][i].vcfg .< cnfg_rw[i])[end] for i in 1:length(cnfg_rw)]
    cnfg_trunc_ts190 = [findall(pp_dat[2][i].vcfg .< cnfg_rw[i])[end] for i in 1:length(cnfg_rw)] ## rwf missing some configs at the end
    cnfg_trunc_ts190_c = [findall(pp_dat_c[2][i].vcfg .< cnfg_rw[i])[end] for i in 1:length(cnfg_rw)]
    truncate_data!(pp_dat[1:2:end], cnfg_trunc_ts001)
    truncate_data!(ap_dat[1:2:end], cnfg_trunc_ts001)
    truncate_data!(pp_dat_c[1:2:end], cnfg_trunc_ts001_c)
    truncate_data!(ap_dat_c[1:2:end], cnfg_trunc_ts001_c)
    truncate_data!(pp_dat[2:2:end], cnfg_trunc_ts190)
    truncate_data!(ap_dat[2:2:end], cnfg_trunc_ts190)
    truncate_data!(pp_dat_c[2:2:end], cnfg_trunc_ts190_c)
    truncate_data!(ap_dat_c[2:2:end], cnfg_trunc_ts190_c)

    if sym_bool[ens.id] == true
        pp = corr_obs_TSM.(pp_dat[1:length(ap_dat_c)], pp_dat_c[1:length(ap_dat_c)], rw=rwf, L=ens.L, info=info, replica_sl=ens.cnfg, nms=sum(ens.cnfg))
        ap = corr_obs_TSM.(ap_dat[1:length(ap_dat_c)], ap_dat_c[1:length(ap_dat_c)], rw=rwf, L=ens.L, info=info, replica_sl=ens.cnfg, nms=sum(ens.cnfg))

        if ens.id == "E250"
            pp_sym = [corr_sym_E250(pp[i], pp[i+1], +1) for i in 1:2:length(ap)]
            ap_sym = [corr_sym_E250(ap[i], ap[i+1], -1) for i in 1:2:length(ap)]
        elseif ens.id == "D450"
            pp_sym = [corr_sym_D450(pp_ts001[i], pp_tsT[i], +1) for i in 1:length(pp_ts001)]
            ap_sym = [corr_sym_D450(ap_ts001[i], ap_tsT[i], -1) for i in 1:length(pp_ts001)]
        else
            pp_sym = [corr_sym(pp[i], pp[i+1], +1) for i in 1:2:length(ap)]
            ap_sym = [corr_sym(ap[i], ap[i+1], -1) for i in 1:2:length(ap)]
        end
    else
        pp_ts001 = corr_obs_TSM.(pp_dat[1:2:length(ap_dat)], pp_dat_c[1:length(ap_dat_c)], rw=rwf, L=ens.L, info=info, replica_sl=ens.cnfg, nms=sum(ens.cnfg))
        ap_ts001 = corr_obs_TSM.(ap_dat[1:2:length(ap_dat)], ap_dat_c[1:length(ap_dat_c)], rw=rwf, L=ens.L, info=info, replica_sl=ens.cnfg, nms=sum(ens.cnfg))
        pp_tsT = corr_obs.(pp_dat[2:2:length(ap_dat)], rw=rwf, L=ens.L, info=info, replica=ens.cnfg, nms=sum(ens.cnfg))
        ap_tsT = corr_obs.(ap_dat[2:2:length(ap_dat)], rw=rwf, L=ens.L, info=info, replica=ens.cnfg, nms=sum(ens.cnfg))

        if ens.id == "E250"
            pp_sym = [corr_sym_E250(pp_ts001[i], pp_tsT[i], +1) for i in 1:length(pp_ts001)]
            ap_sym = [corr_sym_E250(ap_ts001[i], ap_tsT[i], -1) for i in 1:length(pp_ts001)]
        elseif ens.id == "D450"
            pp_sym = [corr_sym_D450(pp_ts001[i], pp_tsT[i], +1) for i in 1:length(pp_ts001)]
            ap_sym = [corr_sym_D450(ap_ts001[i], ap_tsT[i], -1) for i in 1:length(pp_ts001)]
        else
            pp_sym = [corr_sym(pp_ts001[i], pp_tsT[i], +1) for i in 1:length(pp_ts001)]
            ap_sym = [corr_sym(ap_ts001[i], ap_tsT[i], -1) for i in 1:length(pp_ts001)]
        end
    end

    return pp_sym, ap_sym
end

function get_corr_TSM(path::String, ens::EnsInfo, g1::String, g2::String; rw=false, info=false, legacy=false, fs=false)
    path_rw = joinpath(path, ens.id, "rwf_def")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    if length(path_rw) == 0
        path_rw = joinpath(path, ens.id, "rwf")
        global path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    end
    path_sl = joinpath(path, ens.id, "sloppy")
    path_sl = filter(x->occursin(".mesons.dat", x), readdir(path_sl, join=true))
    path_c = joinpath(path, ens.id, "correc")
    path_c = filter(x->occursin(".mesons.dat", x), readdir(path_c, join=true))

    rwf = read_ms1.(path_rw, v=ens.vrw)
    dat_sl = read_mesons([path_sl[i] for i in 1:length(path_sl)], g1, g2, legacy=legacy, id=ens.id)
    dat_c = read_mesons_correction([path_c[i] for i in 1:length(path_c)], g1, g2, legacy=legacy, id=ens.id)

    rw ? corr = [corr_obs_TSM(dat_sl[i], dat_c[i], L=ens.L, rw=rwf, info=info) for i in 1:length(dat_sl)] : corr = [corr_obs_TSM(dat_sl[i], dat_c[i], L=ens.L, info=info) for i in 1:length(dat_sl)]

    if info == false
        return corr
    else
        pp = getindex.(corr,1)
        ppw = getindex.(corr,2)
        w = getindex.(corr,3)
        return pp, ppw, w[1]
    end
end

function get_corr_wil(path::String, ens::EnsInfo, g1::String, g2::String; rw=false, info=false, legacy=false, fs=false)
    path_rw = joinpath(path, ens.id, "rwf_def")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    if length(path_rw) == 0
        path_rw = joinpath(path, ens.id, "rwf")
        global path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    end
    path = joinpath(path, ens.id, "wil")
    path = filter(x->occursin(".mesons.dat", x), readdir(path, join=true))

    if ens.id == "D200"
        #rwf_1 = read_ms1.([path_rw[1]], v=ens.vrw)
        #rwf_2 = read_ms1.([path_rw[2]], v=ens.vrw)
        #rwf = [hcat(rwf_1[1],rwf_2[1])]
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat = read_mesons([path[1]], g1, g2, legacy=legacy, id=ens.id)
        dat_2 = read_mesons([path[2]], g1, g2, legacy=legacy, id=ens.id)
        concat_data!(dat,dat_2)
    else
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat = read_mesons([path[i] for i in 1:length(path)], g1, g2, legacy=legacy, id=ens.id)
        truncate_data!(dat,ens.cnfg)
    end

    rw ? corr = [corr_obs(dat[i], L=ens.L, rw=rwf, info=info, flag_strange=fs) for i in 1:length(dat)] : corr = [corr_obs(dat[i], L=ens.L, info=info) for i in 1:length(dat)]

    if info == false
        return corr
    else
        pp = getindex.(corr,1)
        ppw = getindex.(corr,2)
        w = getindex.(corr,3)
        return pp, ppw, w[1]
    end
end

function get_corr_tm(path::String, ens::EnsInfo, g1::String, g2::String; rw=false, info=false, legacy=false, fs=false)
    path_rw = joinpath(path, ens.id, "rwf_def")
    path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    if length(path_rw) == 0
        path_rw = joinpath(path, ens.id, "rwf")
        global path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
    end
    path = joinpath(path, ens.id, "tm")
    path = filter(x->occursin(".mesons.dat", x), readdir(path, join=true))

    if ens.id == "J303"
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat = read_mesons([path[1]], g1, g2, legacy=legacy, id=ens.id)
        dat_2 = read_mesons([path[2]], g1, g2, legacy=legacy, id=ens.id)
        concat_data!(dat,dat_2)
    elseif ens.id == "D200"
        #rwf_1 = read_ms1.([path_rw[1]], v=ens.vrw)
        #rwf_2 = read_ms1.([path_rw[2]], v=ens.vrw)
        #rwf = [hcat(rwf_2[1],rwf_1[1])]
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat = read_mesons([path[1]], g1, g2, legacy=legacy, id=ens.id)
        dat_2 = read_mesons([path[2]], g1, g2, legacy=legacy, id=ens.id)
        dat_3 = read_mesons([path[3]], g1, g2, legacy=legacy, id=ens.id)
        concat_data!(dat,dat_3)
        concat_data!(dat,dat_2)
    elseif ens.id == "N300"
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat_1 = read_mesons([path[1]], g1, g2, legacy=legacy, id=ens.id)
        dat = read_mesons([path[2]], g1, g2, legacy=legacy, id=ens.id)
        truncate_data!(dat,[199])
        concat_data!(dat,dat_1)
    else
        rwf = read_ms1.(path_rw, v=ens.vrw)
        dat = read_mesons([path[i] for i in 1:length(path)], g1, g2, legacy=legacy, id=ens.id)
        truncate_data!(dat,ens.cnfg)
    end

    rw ? corr = [corr_obs(dat[i], L=ens.L, rw=rwf, info=info) for i in 1:length(dat)] : corr = [corr_obs(dat[i], L=ens.L, info=info) for i in 1:length(dat)]

    if info == false
        return corr
    else
        pp = getindex.(corr,1)
        ppw = getindex.(corr,2)
        w = getindex.(corr,3)
        return pp, ppw, w[1]
    end
end

function read_ens_TSM(path::String, ens::EnsInfo; legacy=false, fs=false)
    pp = get_corr_TSM(path, ens, "G5", "G5", rw=true, info=false, legacy=legacy);
    ap = get_corr_TSM(path, ens, "G5", "G0G5", rw=true, info=false, legacy=legacy);
    idx_wil = findall([pp[i].mu == [.0,.0] for i in 1:length(pp)]);
    idx_tm = findall([pp[i].mu[1] != .0 for i in 1:length(pp)]);

    pp_sym = [corr_sym(pp[i], pp[i+1], +1) for i in idx_wil[1]:2:idx_wil[end]-1];
    ap_sym = [corr_sym(ap[i], ap[i+1], -1) for i in idx_wil[1]:2:idx_wil[end]-1];
    pp_tm_sym = [corr_sym(pp[i], pp[i+div(length(idx_tm),2)], +1) for i in idx_tm[1]:idx_tm[1]+div(length(idx_tm),2)-1];
    ap_tm_sym = [corr_sym(ap[i], ap[i+div(length(idx_tm),2)], -1) for i in idx_tm[1]:idx_tm[1]+div(length(idx_tm),2)-1];

    pp_sym = [pp_sym[1:3]; pp_tm_sym]
    ap_sym = [ap_sym[1:3]; ap_tm_sym]

    return pp_sym, ap_sym
end

function read_ens_wil(path::String, ens::EnsInfo; legacy=false, fs=false)
    pp, ppw, w = get_corr_wil(path, ens, "G5", "G5", rw=true, info=true, legacy=legacy);
    pp_sym = [corr_sym(pp[i], pp[i+1], +1) for i in 1:2:length(pp)-1];
    ap, apw, w = get_corr_wil(path, ens, "G5", "G0G5", rw=true, info=true, legacy=legacy);
    ap_sym = [corr_sym(ap[i], ap[i+1], -1) for i in 1:2:length(ap)-1];

    pp_d1 = get_corr_wil(path, ens, "G5_d1", "G5_d1", rw=true, legacy=legacy);
    pp_d2 = get_corr_wil(path, ens, "G5_d2", "G5_d2", rw=true, legacy=legacy);
    ap_d1 = get_corr_wil(path, ens, "G5_d1", "G0G5_d1", rw=true, legacy=legacy);
    ap_d2 = get_corr_wil(path, ens, "G5_d2", "G0G5_d2", rw=true, legacy=legacy);
    if ens.id == "D200"
        dSdm = get_dSdm(path, ens)
        aux = [hcat(dSdm[3], dSdm[1])]
        aux = [hcat(aux[1], dSdm[2])]
        dSdm = aux
    else
        dSdm = get_dSdm(path, ens)
    end

    pp_val = [[pp_d1[i], pp_d2[i]] for i in 1:length(pp_d1)];
    ap_val = [[ap_d1[i], ap_d2[i]] for i in 1:length(ap_d1)];
    corr = [[pp[i] for i in 1:length(pp)]; [ap[i] for i in 1:length(ap)]];
    corr_val = [[pp_val[i] for i in 1:length(pp)]; [ap_val[i] for i in 1:length(ap)]];
    corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

    return pp_sym, ap_sym, corr, corr_val, corrw, dSdm, w
end

function read_ens_tm_sym(path::String, ens::EnsInfo; legacy=false)
    pp, ppw, w = get_corr_tm(path, ens, "G5", "G5", rw=true, info=true, legacy=legacy);
    pp_sym = [corr_sym(pp[i], pp[i+9], +1) for i in 1:9];
    ap, apw, w = get_corr_tm(path, ens, "G5", "G0G5", rw=true, info=true, legacy=legacy);
    ap_sym = [corr_sym(ap[i], ap[i+9], -1) for i in 1:9];

    dSdm = get_dSdm(path, ens)

    corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

    return pp_sym, ap_sym, corrw, dSdm, w
end

function read_ens_tm(path::String, ens::EnsInfo; legacy=false)
    pp, ppw, w = get_corr_tm(path, ens, "G5", "G5", rw=true, info=true, legacy=legacy);
    pp_sym = [corr_sym(pp[i], pp[i+24], +1) for i in 1:24];
    ap, apw, w = get_corr_tm(path, ens, "G5", "G0G5", rw=true, info=true, legacy=legacy);
    ap_sym = [corr_sym(ap[i], ap[i+24], -1) for i in 1:24];

    if ens.id == "D200"
        dSdm = get_dSdm(path, ens)
        aux = [hcat(dSdm[3], dSdm[1])]
        aux = [hcat(aux[1], dSdm[2])]
        dSdm = aux
    else
        dSdm = get_dSdm(path, ens)
    end

    corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

    return pp_sym, ap_sym, corrw, dSdm, w
end

function read_ens_csv(ens::EnsInfo)
    ix = ensemble_inv[ens.id]
    path_ll = Path_ll[ix]
    path_ls = Path_ls[ix]
    path_ss = Path_ss[ix]
    path_ll_sym = Path_ll_sym[ix]
    path_ls_sym = Path_ls_sym[ix]
    path_ss_sym = Path_ss_sym[ix]
    path2_ll = Path2_ll[ix]
    path2_ls = Path2_ls[ix]
    path2_ss = Path2_ss[ix]
    path2_ll_sym = Path2_ll_sym[ix]
    path2_ls_sym = Path2_ls_sym[ix]
    path2_ss_sym = Path2_ss_sym[ix]

    pp_ll = [csv2Corr(path_ll[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ll)]
        pp_ls = [csv2Corr(path_ls[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ls)]
    pp_ss = [csv2Corr(path_ss[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ss)]
    ap_ll = [csv2Corr(path2_ll[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ll)]
        ap_ls = [csv2Corr(path2_ls[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ls)]
    ap_ss = [csv2Corr(path2_ss[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ss)]

    pp_ll_2 = [csv2Corr(path_ll_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ll_sym)]
        pp_ls_2 = [csv2Corr(path_ls_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ls_sym)]
    pp_ss_2 = [csv2Corr(path_ss_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path_ss_sym)]
    ap_ll_2 = [csv2Corr(path2_ll_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ll_sym)]
        ap_ls_2 = [csv2Corr(path2_ls_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ls_sym)]
        ap_ss_2 = [csv2Corr(path2_ss_sym[i], ens_cnfg[ens.id], trunc=ens.cnfg, id=ens.id) for i in 1:length(path2_ss_sym)]

        pp_ll_sym = [corr_sym(pp_ll[i],pp_ll_2[i],+1) for i in 1:1:length(pp_ll)]
        pp_ls_sym = [corr_sym(pp_ls[i],pp_ls_2[i],+1) for i in 1:1:length(pp_ls)]
        ap_ll_sym = [corr_sym(ap_ll[i],ap_ll_2[i],-1) for i in 1:1:length(ap_ll)]
        ap_ls_sym = [corr_sym(ap_ls[i],ap_ls_2[i],-1) for i in 1:1:length(ap_ls)]
        pp_ss_sym = [corr_sym(pp_ss[i],pp_ss_2[i],+1) for i in 1:1:length(pp_ss)]
        ap_ss_sym = [corr_sym(ap_ss[i],ap_ss_2[i],-1) for i in 1:1:length(ap_ss)]

    i=j=1
    pp_sym = Array{juobs.Corr,1}()
    ap_sym = Array{juobs.Corr,1}()
    while i < length(pp_ll_sym)
        pp_sym = vcat(pp_sym, [pp_ll_sym[i:i+2]; pp_ls_sym[j:j+8]; pp_ss_sym[i:i+2]])
        ap_sym = vcat(ap_sym, [ap_ll_sym[i:i+2]; ap_ls_sym[j:j+8]; ap_ss_sym[i:i+2]])
        i+=3
        j+=9
    end

    return pp_sym, ap_sym, [pp_ll;pp_ls;pp_ss;pp_ll_2;pp_ls_2;pp_ss_2], [ap_ll;ap_ls;ap_ss;ap_ll_2;ap_ls_2;ap_ss_2]
end
=#
