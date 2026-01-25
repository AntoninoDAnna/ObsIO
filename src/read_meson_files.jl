using Statistics, ADerrors

function read_GlobalHeader(path::String)::GlobalHeader
    data =   open(path, "r")
    gh =  zeros(Int32, 4)
    read!(data, gh)
    close(data)
    return GlobalHeader(gh)
end

function __read_GlobalHeader(file)::GlobalHeader
    gh = zeros(Int32,4)
    read!(file, gh)
    return GlobalHeader(gh)
end

function __read_legacy_CorrHeader!(ch,file)::Nothing
    aux_f = zeros(Float64, 4)
    aux_i = zeros(Int32, 4)
    @inbounds for k in eachindex(ch)
        read!(file, aux_f)
        read!(file, aux_i)
        ch[k] = CorrHeader(aux_f, aux_i)
    end
end

function __read_Qsmearing(file)::Smearing{QuarkSmearing.Type}
    qs = read(file,Int32) |> QuarkSmearing.Type
    qs == QuarkSmearing.Local && (return Smearing(qs))
    n = read(file,Int32)
    eps = read(file,Float64)
    return Smearing(qs,n,eps)
end

function __read_Gsmearing(file,qs::Smearing{QuarkSmearing.Type})::Smearing{GluonicSmearing.Type}
    gs = read(file,Int32) |> GluonicSmearing.Type
    gs == GluonicSmearing.Local && (return Smearing(gs))
    if (gs == GluonicSmearing.Quark3DGradientFlow ||
        gs == GluonicSmearing.QuarkGradientFlow)
        return Smearing(gs,qs.niter,qs.eps)
    end
    n = read(file,Int32)
    eps = read(file,Float64)
    return Smearing(gs,n,eps)
end

function __read_CorrHeader!(ch,file)::Nothing
    aux_f = zeros(Float64, 6)
    aux_i = zeros(Int32, 4)
    theta = zeros(Float64, 6)
    @inbounds for k in eachindex(ch)
        read!(file, aux_f)
        read!(file, theta)
        qs1 = __read_Qsmearing(file)
        qs2 = __read_Qsmearing(file)
        gs1 = __read_Gsmearing(file,qs1)
        gs2 = __read_Gsmearing(file,qs2)
        read!(file, aux_i)
        ch[k] =  CorrHeader(aux_f, aux_i, theta, [qs1, qs2, gs1, gs2])
    end
end

function read_CorrHeader(path::String; legacy::Bool=false)
    data = open(path, "r")
    gh = __read_GlobalHeader(data)
    a = Vector{CorrHeader}(undef, gh.ncorr)
    legacy ? __read_legacy_CorrHeader!(a,data) : __read_CorrHeader!(a,data)
    close(data)
    return a
end

function __find_match(ch,g1::Gamma,g2::Gamma)
    function f(x)
        return (x.type[1] == g1 || g1 ==None) &&  (x.type[2] == g2 || g2==None)
    end
    findall(f,ch)
end

function __find_match(ch,gamma::NTuple{2,Gamma}...)
    function f(x)
        return x.type in gamma
    end
    findall(f,ch)
end

function __read_mesons(path,gh,ch,match;
                      nnoise_trunc=false,
                      legacy=false,
                      id,
                      correction::Bool = false)
    ncorr = gh.ncorr
    tvals = gh.tvals
    nnoise = gh.nnoise
    nnoise_trunc = isnothing(nnoise_trunc) ? nnoise : min(nnoise, nnoise_trunc)
    ncfg = let
        fsize = filesize(path)
        datsize =  4 + sum(c.dsize for c in ch) * tvals * nnoise
        div(fsize - gh.hsize - sum(c.hsize for c in ch), datsize)
    end
    N = correction ? div(length(match),2) : length(match)
    data_re = zeros(N, ncfg, tvals)
    data_im = zeros(N, ncfg, tvals)
    vcfg  = Array{Int32}(undef, ncfg)
    tmp   = Array{Float64}(undef, 2*tvals*nnoise)
    tmp_v = view(tmp,1:tvals*nnoise)
    open(path,"r") do data
        seek(data, gh.hsize + sum(c.hsize for c in ch))
        @inbounds for icfg = 1:ncfg
            vcfg[icfg] = read(data, Int32)
            c,sgn=1,+1 ## if correction, at div(ncorr,2) it changes to -
            @inbounds for k = 1:ncorr
                if !(k in match)
                    seek(data, position(data)  + ch[k].dsize*tvals*nnoise)
                    continue
                end
                if ch[k].is_real
                    read!(data, tmp_v)
                    tmp2 = reshape(tmp_v, (nnoise, tvals))
                    tmp2 = mean(tmp2[1:nnoise_trunc, :], dims=1)
                    data_re[c, icfg, :] .+= sgn * tmp2[1, :]
                else
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
    @inbounds for c in 1:N
        cm = match[c]
        res[c] = CorrData(ch[cm],vcfg,data_re[c, :, :],data_im[c, :, :],id)
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
                     id = nothing,
                     legacy = false,
                     nnoise_trunc = nothing)
    id = get_id(path,id)
    gh = read_GlobalHeader(path)
    c_header = read_CorrHeader(path, legacy=legacy)
    corr_match = __find_match(c_header,g1,g2)
    return  __read_mesons(path,g_header,c_header,corr_match,
                         nnoise_trunc = nnoise_trunc, legacy=legacy,id=id)
end

function read_mesons(path::String,
                     gamma::NTuple{2,Gamma}...;
                     id = nothing,
                     legacy = false,
                     nnoise_trunc = nothing)
    id = get_id(path,id)
    g_header = read_GlobalHeader(path)
    c_header = read_CorrHeader(path, legacy=legacy)
    corr_match = __find_match(c_header,gamma...)
    return __read_mesons(path,g_header,c_header,corr_match,
                         nnoise_trunc = nnoise_trunc, legacy=legacy,id=id)
end

function read_mesons(path::AbstractVector{String}, p...;k... )
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
                                id = nothing,
                                legacy = false,
                                nnoise_trunc = nothing)
    id = get_id(path,id)
    g_header = read_GlobalHeader(path)
    c_header = read_CorrHeader(path, legacy=legacy)
    match = __find_match(c_header,g1,g2)
    return __read_mesons(path,g_header,c_header,match,id=id,nnoise_trunc=nnoise_trunc,
                        legacy=legacy,correction=true)
end

function read_mesons_correction(path::String,
                                gamma::NTuple{2,Gamma}...;
                                id = nothing,
                                legacy = false,
                                nnoise_trunc = nothing)
    id = get_id(path,id)
    g_header = read_GlobalHeader(path)
    c_header = read_CorrHeader(path, legacy=legacy)
    match = __find_match(c_header,gamma...)
    return __read_mesons(path,g_header,c_header,match,id=id,nnoise_trunc=nnoise_trunc,
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

function read_mesons_by_conf(path::Vector{String},gamma...;
                            nconf = length(path),
                            id=nothing, legacy=false, nnoise_trunc=nothing)
    id = get_id(path[1],id)
    gh = read_GlobalHeader(path[1])
    ch = read_CorrHeader(path[1])
    match = __find_match(ch,gamma...)
    Ncorr = length(match)
    res = [CorrData(ch[i],nconf,gh.tvals,id) for i in eachindex(match)]
    for p in path
        _r = __read_mesons(p,gh,ch,match,id=id,legacy=legacy,
                                nnoise_trunc=nnoise_trunc)
        for i in eachindex(_r)
            for cdx in eachindex(_r[i].vcfg)
                res[i].re_data[_r[i].vcfg[cdx],:] .= _r[i].re_data[cdx,:]
                res[i].im_data[_r[i].vcfg[cdx],:] .= _r[i].im_data[cdx,:]
            end
        end
    end
    return res
end

function apply_rw(data::AbstractArray{Float64}, W::AbstractMatrix{Float64}, vcfg = nothing; id = nothing, fs = false)
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

function apply_rw(data::AbstractVector{<:AbstractArray{Float64}},
                  W::AbstractVector{<:AbstractMatrix{Float64}}, vcfg = nothing; id = nothing, fs = false)
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
                  real = true,
                  rw::Union{Array{Float64, 2}, Nothing}=nothing,
                  L = 1, info = false,
                  idm = nothing,
                  nms = Int64(maximum(cdata.vcfg)),
                  flag_strange = false)
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
    corr = __update__(corr,obs=obs)
    if info
        return !isnothing(rw) ?  (corr,obs) : (corr,ow,W_obs)
    else
        return corr
    end
end

function corr_obs(cdata::AbstractVector{CorrData}, corr;
                  real = true,
                  replica = nothing,
                  rw::Union{Vector{Array{Float64, 2}}, Nothing}=nothing,
                  L = 1, info = false,
                  idm = nothing,
                  nms = 0,
                  flag_strange = false)
    nrep = length(cdata)
    id = let
        ids = getfield.(cdata,:id)
        !all(ids .== ids[1]) && error("[corr_obs] IDs are not equal")
        ids[1]
    end
    vcfg = getfield.(cdata,:vcfg)
    replica = isnothing(replica) ? Int64.(maximum.(vcfg)) : replica
    nms==0 && (nms=sum(replica))
    if isnothing(idm)
        a = vcfg[1]
        for i in 2:nrep
            a = [a; a[end] .+ vcfg[i]]
        end
        idm = Int64.(a)
    end
    data = let s = real ? :re_data : :im_data
        getfield.(cdata,s)./L^3
    end
    nt = size(data[1])[2]
    if isnothing(rw)
        obs = vcat(data...) |>
            x->[uwreal(x[:, x0], id, replica, idm, nms) for x0 = 1:nt]
    else
        data_r, W = apply_rw(data, rw, vcfg, id=id, fs=flag_strange)
        ow =vcat(data_r...) |>
            x-> [uwreal(x[:, x0],id, replica, idm, nms) for x0 = 1:nt]
        W_obs =vcat(W...) |>
            x->uwreal(x, id,replica, idm, nms)
        obs = [ow[x0] / W_obs for x0 = 1:nt]
    end
    corr = __update__(corr,obs=obs)
    if info
        return !isnothing(rw) ?  (corr,obs) : (corr,ow,W_obs)
    else
        return corr
    end
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

function concat_data!(data1::Vector{CorrData}, data2::Vector{CorrData})
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

function concat_data!(data1::Vector{Vector{CorrData}}, data2::Vector{Vector{CorrData}})
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
