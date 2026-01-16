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
    seq_prop = -1
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
        elseif contains("mom",line)
            pF = extract_all(r"[0-9]+",line,Int64)
            println(pF)
        elseif contains("seq_x0",line)
            seq_x0 = extract_first(r"[0-9]+",line,Int64)
        elseif  !isnothing(match(r"(?<=\[)([^\]]+)",line))
            error("Overflowed Propagator section")
        end
    end
    prop = Propagator(kappa,mu,theta,pF,-1,-1)
    src  = Point(seq_type,seq_x0,QuarkSmearing.Local,GluonicSmearing.Local)
    snk  = Point(None,missing, QuarkSmearing.Local,GluonicSmearing.Local)
    return prop, src, snk
end

function read_correlator(file,props,points)
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
            type = extract_all(r"(G[0-5]){1,2}",line, Gamma)
        elseif contains(line,"gsmear")
            gsmear = extract_all(r"[0-9]+",line,Int32) |>
                x->GluonicSmearing.Type.(x)
        elseif contains(line,"qsmear")
            qsmear = extract_all(r"[0-9]+",line,Int32) |>
                x->QuarkSmearing.Type.(x)
        elseif contains(line,"x0")
            x0 = extract_first(r"[0-9]+",line,Int64)
        end
    end
    Pr1 = update(props[iprop[1]], src=1)
    src = points[2iprop[1]-1]
    snk = points[2iprop[1]]

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
    points=nothing;
    corrs =nothing;
    while !eof(file)
        head = match(r"(?<=\[)([^\]]+)",readline(file))
        isnothing(head) && continue;
        if contains(head.match,"Measurements")
            nprop,ncorr = read_measurements(file)
            props = Vector{Propagator}(undef, nprop)
            points = Vector{Point}(undef,2nprop)
            corrs = Vector{Corr}(undef,ncorr)
        elseif contains(head.match,"Propagator")
            idx = match(r"[0-9]+",head.match) |> x->parse(Int64,x.match) + 1
            props[idx], points[2idx-1], points[2idx] = read_propagator(file)
        elseif contains(head.match,"Correlator")
            idx = match(r"[0-9]+",head.match)|> x->parse(Int64,x.match) + 1
            read_correlator(file,props,points)
        end
    end
    close(file)
   # return props, points
end
