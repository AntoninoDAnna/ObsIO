function extract_first(regex,line,type)
    m = match(regex,line)
    isnothing(m) && error("Cannot match regex to line \"",line,"\"")
    return parse(type,m.match)
end

function extract_all(regex,line,type)
    m = eachmatch(regex,line)
    isnothing(m) && error("Cannot match regex to line \"",line,"\"")
    return tuple([parse(type,x.match) for x in m]...)
end

function read_propagator(file)
    kappa, mu = 0.0, 0.0
    pF = ntuple(x->0.0,Val(4))
    theta = ntuple(x->0.0,Val(3))
    seq_prop = false
    seq_type = None
    seq_x0 = -1
    flag::UInt8 = 0x0
    while true
        line = readline(file)
        (line == "") && break
        aux = findfirst("#", line)
        if !isnothing(aux)
            line = line[1:aux[1]-1]
        end
        line=="" && continue;
        if contains(line,"kappa")
            kappa = extract_first(r"[+-]?([0-9]*[.])?[0-9]+",line,Float64)
            flag |= 0x1
        elseif contains(line,"mus")
            mu =  extract_first(r"[+-]?([0-9]*[.])?[0-9]+",line,Float64)
            flag |= 0x2
        elseif contains(line,"theta")
            theta = extract_all(r"[+-]?([0-9]*[.])?[0-9]+",line,Float64)
            flag |= 0x04
        elseif contains(line,"seq_prop")
            seq_prop=extract_first(r"[0-9]+",line,Int64)
        elseif contains(line, "seq_type")
            seq_type = extract_first(r"(G[0-5]){1,2}",line,Gamma)
        elseif contains(line,"mom")
            pF = extract_all(r"[0-9]+",line,Int64)
        elseif contains(line,"seq_x0")
            seq_x0 = extract_first(r"(?<=\s)[0-9]+",line,Int64)
        elseif  !isnothing(match(r"(?<=\[)([^\]]+)",line))
            error("Overflowed Propagator section")
        end
    end
    src  = Point(seq_type,seq_x0,QuarkSmearing.Local,GluonicSmearing.Local)
    snk  = Point(None,missing, QuarkSmearing.Local,GluonicSmearing.Local)
    prop = Propagator(kappa,mu,theta,pF,src,snk,seq_prop)
    return prop
end

function resolve_seq_prop(p1::Propagator,p2::Propagator,
                          props::Vector{Propagator})
    res_props = (p1,p2)
    if isa(p1.seq_prop, Int64)
        idx = p1.seq_prop +1
        update(p1,seq_prop = true)
        seq_props = resolve_seq_prop(props[idx],p2,props)
        res_props = (p1,seq_prop...)
    elseif isa(p2.seq_prop,Int64)
        idx = p2.seq_prop +1
        update(p2,seq_prop = true)
        seq_prop = resolve_seq_prop(p1,props[idx],props)
        res_props = (seq_prop...,p2)
    end
        return res_props
end

function update_propagators_tuple(props::NTuple{2,Propagator},
                                  x0, types,qsmear,gsmear)
    p1,p2 = props
    p1 = update(p1,src = update(p1.src,x0=x0, qsmearing = qsmear[1],
                                gsmearing = gsmear[1],gamma = types[1]),
                snk = update(p1.snk,x0=missing, qsmearing = qsmear[2],
                             gsmearing = gsmear[2],gamma = types[2]))
    p2 = update(p2,src = p1.snk)
    p2 = update(p2,snk = p1.src)
    return (p1,p2)
end

function update_propagators_tuple(props::NTuple{3,Propagator},
                                  x0, types,qsmear,gsmear)
    p1,p2,p3 = props # p3 is the sequential propagator

    p1 = update(p1,src = update(p1.src,x0=x0,qsmearing=qsmear[1],
                                gsmearing = gsmear[1],gamma = types[1]),
                snk = update(p1.snk,x0=missing,qsmearing=qsmear[2],
                             gsmearing = gsmear[2],gamma = types[2]))
    p2 = update(p2,src = p1.src)
    p2 = update(p2,snk = p3.src)
    p3 = update(p3,snk = p1.snk)
    return (p1,p2,p3)
end


function read_correlator(file,props)
    iprop = (0.,0.)
    type = nothing
    gsmear = (0,0)
    qsmear = (0,0)
    x0 = -1
    while true
        line = readline(file)
        (line == "") && break
        aux = findfirst("#", line)
        if !isnothing(aux)
            line = line[1:aux[1]-1]
        end
        line=="" && continue;
        if contains(line,"iprop")
            iprop = extract_all(r"[0-9]+",line,Int64) |>
                x->(x[1]+1,x[2]+1)
        elseif contains(line,"type")
            type = extract_all(r"(G[0-5]){1,2}|(1)",line, Gamma)
        elseif contains(line,"gsmear")
            gsmear = extract_all(r"[0-9]+",line,Int32) |>
                x->GluonicSmearing.Type.(x)
        elseif contains(line,"qsmear")
            qsmear = extract_all(r"[0-9]+",line,Int32) |>
                x->QuarkSmearing.Type.(x)
        elseif contains(line,"x0")
            x0 = extract_first(r"(?<=\s)[0-9]+",line,Int64)
        end
    end
    props = resolve_seq_prop(props[iprop[1]],props[iprop[2]],props)
    props = update_propagators_tuple(props, x0,type,qsmear,gsmear)
    return Corr([],props)
end
function read_measurements(file)
    nprop,ncorr = 0,0
    flag::UInt8 = 0x0
    while true
        line = readline(file)
        if contains(line,"nprop")
            nprop = match(r"[0-9]+",line) |> x->parse(Int64, x.match)
            flag |= 0x1
        elseif contains(line,"ncorr")
            ncorr = match(r"[0-9]+",line) |> x->parse(Int64, x.match)
            flag |= 0x2
        end
        ((flag & 0x03) == 0x03) && break;
        if !isnothing(match(r"(?<=\[)([^\]]+)",line))
            error("In input file, section [Measurements], nprop or ncorr are missing")
        end
    end
    return nprop,ncorr
end

function read_input_file(path)
    file = open(path,"r")
    props =nothing;
    corrs =nothing;
    while !eof(file)
        head = match(r"(?<=\[)([^\]]+)",readline(file))
        isnothing(head) && continue;
        if contains(head.match,"Measurements")
            nprop,ncorr = read_measurements(file)
            props = Vector{Propagator}(undef, nprop)
            corrs = Vector{Corr}(undef,ncorr)
        elseif contains(head.match,"Propagator")
            idx = match(r"[0-9]+",head.match) |> x->parse(Int64,x.match) + 1
            props[idx] = read_propagator(file)
        elseif contains(head.match,"Correlator")
            idx = match(r"[0-9]+",head.match)|> x->parse(Int64,x.match) + 1
            corrs[idx] = read_correlator(file,props)
        end
    end
    close(file)
    return corrs
end

function read_input_file(path,gsrc::Gamma, gsnk::Gamma)
    corrs = read_input_file(path)
    gsrc != None && filter!(x->x.points[1].gamma ==gsrc, corrs)
    gsrc != None && filter!(x->x.points[2].gamma==gsnk,corrs)
    return corrs
end

"No Gamma::None allowed in here"
function read_input_file(path,gamma::NTuple{2,Gamma}...)
    corrs = read_input_file(path)
    filter!(x->(x.points[1].gamma,x.points[2].gamma) in gamma,corrs)
    return corrs
end
