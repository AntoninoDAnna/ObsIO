using ALPHAio, ADerrors, BDIO

@doc raw"""
     read_corr_ALPHA(path)::juobs.Corr

Reads the correlator saved in `path` and returns a `juobs.Corr`.
If a `extra` field is present with the correlator parameters, it uses them
to create the `juobs.Corr`, otherwise it set all the parameters to default
values
See also [`DEF`](@ref)
"""
function read_corr_ALPHA(path)::juobs.Corr
    file = BDIO_open(path,"r")
    BDIO_seek!(file)
    info = ALPHAdobs_read_parameters(file)
    obs = Vector{ADerrors.uwreal}()
    for _ = 1:info["size"][1]
        push!(obs,ALPHAdobs_read_next(file))
    end
    BDIO_close!(file)
    extra = get(info,"extra",nothing)
    if isnothing(extra)
        return juobs.Corr(obs);
    end
    kappa = get(extra,"kappa",DEF.mass)
    mu = get(extra,"mu",DEF,mass)
    gamma = get(extra,"gamma",DEF.gamma)
    y0 = get(extra,"gamma",DEF.y0)
    theta1 = get(extra,"theta1",DEF.theta)
    theta2 = get(extra,"theta2",DEF.theta)

    return juobs.Corr(obs,kappa,mu,gamma,y0,theta1,theta2)
end



