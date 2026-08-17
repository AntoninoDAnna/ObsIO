using BDIO, ALPHAio
@doc raw"""
     open_data_file(path,mode; override=false, comment="")

Open a bdio file in `mode` and return the handle. If a file is created, a path to the file is also made

# Modes
- `"r"`: read mode. Throws an error if the file does not exist
- `"w"`: write mode. if the kwargs `override` is set to `true`, it warns the user and delete any existing file. otherwise throws an error if it exists. `comment`, is passed on to `ALPHAdobs_create`
- `"d"`: delete mode. If open the file and override any existing file without warning the user.  `comment`, is passed on to `ALPHAdobs_create`
 """
function open_data_file(path,mode; comment="",override=false)
    mode == "r" && (return BDIO_open(path,mode))
    if isfile(path)
        if mode == "w"  && override
            @warn "overriding $(basename(path))"
            Base.rm(path)
        elseif mode =="d"
            Base.rm(path,force=true)
        else
            error("$(basename(path)) already exists")
        end
    end
    mkpath(dirname(path))
    return ALPHAdobs_create(path,IOBuffer(comment))
end

@doc """
    close_data_file(fb)

close the data file
"""
close_data_file(fb) = BDIO_close!(fb)

@doc """
     save_data(a,path::String,comment::AbstractString=""; override::Bool = false, extra=nothing)

saves the data `a` in a BDIO file at `path` following the ALPHA conventions. If `override` is given, then
"""
function save_data(a, path::String, comment::AbstractString=""; override::Bool = false, extra=nothing)
    file = open_data_file(path,"w",comment=comment,override=override)
    save_data(a,file,extra=extra)
    ALPHAdobs_close(file)
end

function save_data(a::Dict{String, <:AbstractArray{uwreal}}, fb; extra=nothing)
    kks = collect(keys(a))
    sz =  size(a[kks[1]])
    ALPHAio.ALPHAdobs_vdict_parameters(fb,1, kks, sz; extra=extra)
    for k in values(a)
        ALPHAdobs_add(fb, collect(k))
    end
end

function save_data(a::AbstractArray{<:AbstractArray{uwreal}}, fb; extra=nothing)
    save_data(collect(a[1]),fb,extra=extra)
    for v in a[2:end]
        save_data(collect(v),fb,extra=nothing)
    end
end

save_data(a, fb; extra=nothing) = ALPHAdobs_write(fb,a,extra=extra)

save_data(a;path::String,comment::AbstractString="", override::Bool=false,extra = nothing) = save_data(a,path,comment,override=override, extra=extra)


function get_extras(path)
    file = open_data_file(path,"r")
    extra = Vector{Dict{String,Any}}()
    while ALPHAdobs_next_p(file)
        info = ALPHAdobs_read_parameters(file)
        push!(extra,get(info,"extra",nothing))
    end
    return length(extra[1])==1 ? extra[1] : extra
end

function read_data(path::AbstractString; get_extra = false,
                   read = ADerrors.read_uwreal)
    file =open_data_file(path,"r")
    res = Any[]
    while ALPHAdobs_next_p(file)
        info =  ALPHAdobs_read_parameters(file)

        type = info["type"]
        nobs = info["nobs"]
        dim  = info["dimensions"]
        size = dim ==0 ?  0 : info["size"]

        if type == "obs"
            if dim == 0
                data = uwreal[]
                for i in 1:nobs
                    push!(data,ObsIO_read_next(file,read=read))
                end
            else
                data = [ones(uwreal,size...) for _ in 1:nobs]
                for n in 1:nobs, c in CartesianIndices(data[n])
                    data[n][c] = ObsIO_read_next(file,read=read)
                end
            end
        elseif type == "dict"
            keys = info["keys"]
            if dim ==0
                data = [Dict{String,uwreal}(k=>ObsIO_read_next(file,read=read)  for k in keys) for _ in 1:nobs]
            else
                data = [Dict{String,Array{uwreal,dim}}() for _ in 1:nobs]
                for n in 1:nobs, k in keys
                    data[n][k] = ones(uwreal,size...)
                    for c in CartesianIndices(data[n][k])
                        data[n][k][c] = ObsIO_read_next(file,read=read)
                    end
                end
            end
        end

        if get_extra
            nobs ==1 ? push!(res,(data[1],get(info,"extra",nothing))) : push!(res,(data, get(info,"extra",nothing)))
        else
            nobs ==1 ? push!(res,data[1]) : push!(res,data)
        end
    end

    return length(res) ==1 ?  res[1] : res
end
