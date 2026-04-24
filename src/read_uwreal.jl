import ADerrors, BDIO

"""
        ObsIO_read_uwreal(fb,ws::ADerrors.wspace, mapsids::Dict{Int64, Int64},mapstrids::Dict{String,String},mapvrep::Dict{String,Vector{Int64}})

Read a uwreal object from a bdio record. the file must be in the correct record to read the uwreal. The uwreal is added the `ws` workspace. `mapsid` maps the numeric id to a new numeric id. the new id must no exist already in `ws`.
`mapstrids` maps the string id to a new string id. the new id must not already exists in `ws`. `mapvrep` allows to correct the montecarlo length of a given id. It accept a vector with the new montecarlo lengths, one for `Int64` for each replica. If the new length is larger that what is stored in the bdio, the the montecarlo lengh is extend addign gaps at the end of the montecarlo chain and the delta are rescaled accordingly. If the new length is shorter, a warning is printed and the last configurations are thrown to match the new length.


## example:

temp = uwreal([rand(1000),rand(500)],"temp") # numeric id 123

save_uwrea("temp.bdio",temp);

fb = BDIO_open("temp.bdio","r")
BDIO_seek!(fb);

ObsIO_read_uwreal(fb); #same ad read

##  ObsIO_read_uwreal(fb, Dict(123=>459)); #  read the uwreal and the id 123 get updated to 459
##  ObsIO_read_uwreal(fb, Dict("temp"=>"temp1")) # the id  "temp" get updated to "temp1"
##  ObsIO_read_uwreal(fb, Dict("temp"=>[1100,400])) # now the replica 1 ha 100 gaps at the end, while the last 100 measurement of replica 2 are thrown
"""
function ObsIO_read_uwreal(fb,ws::ADerrors.wspace,mapsids::Dict{Int64,Int64},
                        mapstrids::Dict{String,String},
                        mapvrep::Dict{String,Vector{Int64}})
    ## variable to temporary read data:
    dfoo = zeros(Float64,1)
    ifoo = zeros(Int32,1)
    mean::Float64, nid::Int32 = let
        BDIO.BDIO_read(fb,dfoo)
        BDIO.BDIO_read(fb,ifoo)
        dfoo[1], ifoo[1]
    end
    p = [false for i in 1:ws.nob+nid]
    p[ws.nob+1:end] .= true
    d = zeros(Float64,ws.nob+nid)
    d[ws.nob+1:end] .=1.0
    nds = zeros(Int32,nid)
    BDIO.BDIO_read(fb,nds)
    nrep = zeros(Int32, nid)
    BDIO.BDIO_read(fb, nrep)
    ivrep = zeros(Int32, sum(nrep))
    BDIO.BDIO_read(fb, ivrep)
    ids  = zeros(Int32, nid)
    BDIO.BDIO_read(fb, ids)
    dfl= zeros(Float64,nid)
    itmp = zeros(Int32,nid)
    BDIO.BDIO_read(fb,itmp)
    for i in 1:2
      BDIO.BDIO_read(fb,dfl)
    end
    dfl =Vector{Vector{Float64}}(undef,nid)
    is =1
    for i in 1:nid
       ie = is + nrep[i] - 1
        if (sum(ivrep[is:ie]) != nds[i])
            throw("Replica sum does not match number of measurements")
        end
        dfl[i] = zeros(Float64, nds[i])
        BDIO.BDIO_read(fb, dfl[i])
        is = ie + 1
    end
    ids_obs = Vector{Int64}(undef,nid)
    id2str = Dict{Int64,String}();
    if BDIO.BDIO_eor(fb)
        is = 1
        for i in 1:nid
            ie = is + nrep[i] - 1
            ids_obs[i] = convert(Int64, ids[i]) |> x->get(mapids,x,x)
            id2str[ids[i]] = ADerrors.get_name_from_id(ids_obs[i], ws)
        end
    else
        name = BDIO.BDIO_read_str(fb) #Obs name, aderros throw it away
        is = 1
        for i in 1:nid
            ie = is + nrep[i] - 1
            BDIO.BDIO_read(fb, ifoo)
            str = BDIO.BDIO_read_str(fb)
            id2str[ids[i]] = str ## keep track of the id *inside* the record
            str = get(mapstrids,str,str)
            # Add the string id we want in the database
            ids_obs[i] = ADerrors.get_id_from_name(str, ws)

        end
        is=1
    end

    # at this point we have a complete map of the ids and string ids *inside*
    # the bdio record, so we can reorganize the data and get update the
    # varius ids as we wish, keeping in mind that the string id rules.
    is = 1
    for i in 1:nid
        ie = is +nrep[i]-1
        str = id2str[ids[i]]
        newvrep = get(mapvrep,str,convert(Vector{Int64},ivrep[is:ie]))
        ndfl = zeros(sum(newvrep))
        isn = 1
        iso = 1
        for irep in 1:nrep[i]
            if newvrep[irep] < ivrep[is+irep-1]
                @warn "ID: $(ids[i]), repl $irep. you're throwing data away"
                ien = isn + newvrep[irep]-1
                ieo = iso + newvrep[irep]-1
                ndfl[isn:ien] .= dfl[i][iso:ieo]
                isn = ien+1
                iso = iso + ivrep[is+irep-1]
            else
                ien = isn+ivrep[is+irep-1]-1
                ieo = iso+ivrep[is+irep-1]-1
                ndfl[isn:ien] .= dfl[i][iso:ieo]
                ndfl .= length(ndfl)/length(dfl[i]).*ndfl
                isn = isn + newvrep[irep]
                iso = ieo+1
            end
        end
        # add_DB, with do_maps=false, add the data to the workspace without mapping it to any id, hence, this id is not stored at this point.
        ADerrors.add_DB(ndfl, ids_obs[i], newvrep, ws, true)
        is = ie + 1
    end
    # Here not to use ids[i]
    if BDIO.BDIO_eor(fb)
        is = 1
        for i in 1:nid
            ie = is + nrep[i] - 1
            v = Vector{String}(undef, nrep[i])
            idc = Vector{Int32}(undef, nds[i])
            str = ADerrors.get_name_from_id(convert(Int64, ids[i]), ws)
            for j in 1:nrep[i]
                v[j] = str*"_r"*string(j-1)
            end
            iof = 0
            for j in 1:nrep[i]
                for k in 1:ivrep[is+j-1]
                    idc[k+iof] = convert(Int32, k)
                end
                iof = iof + ivrep[is+j-1]
            end
            ADerrors.add_repnames(convert(Int64, ids_obs[i]), ws, v, convert(Vector{Int64}, idc))
            is = ie + 1
        end
    else
        for i in 1:nid
            v = Vector{String}(undef, nrep[i])
            idc = Vector{Int32}(undef, nds[i])
            ifoo = zeros(Int32,1)
            BDIO.BDIO_read(fb, ifoo)
            for j in 1:nrep[i]
                v[j] = BDIO.BDIO_read_str(fb)
            end
            str = ADerrors.get_name_from_id(ids_obs[i],ws)
            if ids_obs[i]!=ifoo[1]  #if true, ids is changing and repid need to reflect that
                for j in 1:nrep[i]
                    v[j] = str*"_r"*string(j-1)
                end
            end
            BDIO.BDIO_read(fb, idc)
            if (length(idc) == 1)
                ADerrors.add_repnames(convert(Int64, ids_obs[i]), ws, v, [1,2])
            else
                ADerrors.add_repnames(convert(Int64, ids_obs[i]), ws, v, convert(Vector{Int64}, idc))
            end
        end
    end
    return ADerrors.uwreal(mean, p, d)
end



ObsIO_read_uwreal(file) = ObsIO_read_uwreal(file,
                                      ADerrors.wsg,
                                      Dict{Int64,Int64}(),
                                      Dict{String,String}(),
                                      Dict{String,Vector{Int64}}())

ObsIO_read_uwreal(file,m::Dict{String,String}) =
    ObsIO_read_uwreal(file,
                   ADerrors.wsg,
                   Dict{Int64,Int64}(),
                   m,
                   Dict{String,Vector{Int64}}())

ObsIO_read_uwreal(file,m::Dict{String,Vector{Int64}}) =
    ObsIO_read_uwreal(file,
                   ADerrors.wsg,
                   Dict{Int64,Int64}(),
                   Dict{String,String}(),
                   m)

ObsIO_read_uwreal(file,
               mapstrids::Dict{String,String},
               vrep::Dict{String,Vector{Int64}}) =
                   ObsIO_read_uwreal(file,
                                  ADerrors.wsg,
                                  Dict{Int64,Int64}(),
                                  mapstrids,
                                  vrep)


function ObsIO_read_next(fb; size=nothing, keys=nothing, read = ADerrors.read_uwreal)

    if size == nothing && keys == nothing
        ALPHAdobs_next_dobs(fb)
        return read(fb)
    elseif size != nothing && keys == nothing
        aw = Array{ADerrors.uwreal}(undef, size...)
        for i in 1:length(aw)
            ALPHAdobs_next_dobs(fb)
            aw[i] = read(fb)
        end
        return aw
    elseif size == nothing && keys != nothing
        aw = Dict{String, ADerrors.uwreal}()
        for i in 1:length(keys)
            ALPHAdobs_next_dobs(fb)
            aw[keys[i]] = read(fb)
        end
        return aw
    elseif size != nothing && keys != nothing
        aw = Dict{String, Array{ADerrors.uwreal}}()
        for i in 1:length(keys)
            x = Array{ADerrors.uwreal}(undef, size...)
            for j in 1:length(x)
                ALPHAdobs_next_dobs(fb)
                x[j] = read(fb)
            end
            aw[keys[i]] = x
        end
        return aw
    end

    return nothing
end
