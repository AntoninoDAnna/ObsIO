using EnumClasses

@doc raw"""
     @enum Gamma G0 G1 G2 G3 None G5 Id G0G1 G0G2 G0G3 G0G5 G1G2 G1G3 G1G5 G2G3 G2G5 G3G5

Enum list of Dirac's Gamma matrices. `None` take the place of `G4` for backward compatability.

"""
@enum Gamma G0 G1 G2 G3 None G5 Id G0G1 G0G2 G0G3 G0G5 G1G2 G1G3 G1G5 G2G3 G2G5 G3G5

__to_gamma__ = Dict(string(g) => g for g in instances(Gamma))
__to_gamma__["1"] = Id

import Base:parse
parse(T::Type{Gamma},s::AbstractString) = __to_gamma__[s]


@enumclass Noise::Int32 begin
    Z2
    GAUSS
    U1
    Z4
end

"""
    @enumclass QuarkSmearing::Int32  Local Wuppertal GradientFlow3D GradientFlow None

List of Quark Smearing supported
"""

@enumclass QuarkSmearing::Int32 begin
    Local
    Wuppertal
    GradientFlow3D
    GradientFlow
    None
end

__to_qs__ = Dict(string(qs)=>qs for qs in instances(QuarkSmearing.Type))

parse(T::Type{QuarkSmearing.Type},s::AbstractString) = __to_qs__[s]

"""
        @enumclass GluonicSmearing::Int32 Local APE WilsonFlow3D Quark3DGradientFlow QuarkGradientFlow None

List of Gluonic Smearing supported
"""

@enumclass GluonicSmearing::Int32 begin
    Local
    APE
    WilsonFlow3D
    Quark3DGradientFlow
    QuarkGradientFlow
    None
end

__to_gs__ = Dict(string(gs)=>gs for gs in instances(GluonicSmearing.Type))

parse(T::Type{GluonicSmearing.Type},s::AbstractString) = __to_gs__[s]
