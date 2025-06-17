module HydrogenAdsorption

using Reexport

include("coefficient_matrix.jl")
export coefficient_matrix

include("adsorption_systems.jl")
export IsothermParameters, MDAParameters, DAParameters, MaterialProperties, GeometricParameters, OperationalParameters, AdsorptionParameters, MDA_adsorption!
end
