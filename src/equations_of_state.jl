"""
Non differentiated ideal gas equation.

# Inputs:
- `T`: Temperature in the tank / Vector of temperatures at each radial node / K
- `R`: Ideal gas constant / J/(mol·K)
- `M`: Molar mass of the gas / kg/mol
- `ρ`: Density of the gas / kg/m³ (optional, if provided, pressure will be calculated)
- `P`: Pressure of the gas / Pa (optional, if provided, density will be calculated)

# Outputs:
- `ρ`: Density of the gas / kg/m³ if `P` is provided
- `P`: Pressure of the gas / Pa if `ρ` is provided
"""
function ideal_gas_equation(T, R::Float64, M::Float64; ρ=nothing, P=nothing)
    if isnothing(ρ) && isnothing(P)
        error("Either density (ρ) or pressure (P) must be provided.")
    end

    if !isnothing(ρ) && !isnothing(P)
        error("Only one of density (ρ) or pressure (P) should be provided.")
    end
    
    if !isnothing(ρ)
        # Calculate pressure from density
        P = ρ .* R .* T ./ M
        return P
    else
        # Calculate density from pressure
        ρ = P .* M ./ (R .* T)
        return ρ
    end
end

#function ideal_gas_equation(T::Vector{Float64}, R::Float64, M::Float64, ρ::Float64, str::String)
    
#end