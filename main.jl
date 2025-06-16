using Revise
using HydrogenAdsorption

# parameters
# Tank parameters
V = 2.4946e-3 # Volume of tank / m³
L = 0.4 # Length of the tank / m

# Inputs needed for the coefficient_matrix function
R_T = sqrt(V / (pi * L)) # Radius of the tank / m
dr = 0.00025 / 2 # Radial step size / m
k_eff = 0.4304 # Effective thermal conductivity / W/(m·K)
U = 36 # Heat transfer coefficient / W/(m²·K)
T₀ = 281 # Initial temperature of the tank / K
T_air = 281

n_r, r_span, A, Tᵢ, b = HydrogenAdsorption.coefficient_matrix(R_T, dr, k_eff, U, T₀, T_air)

display(A)
display(Tᵢ)
println("n_r: ", n_r)
println("b: ", b)
println("r_span: ", r_span)

# Isotherm parameters
# Modified Dubinin-Astakov Isotherm parameters
α = 3080.0 # Enthalpic Factor / J/mol
β = 18.9 # Entropic Factor / J/mol K
m = 2.0 # Exponential factor
p₀ = 1470e6 # Saturation pressure / Pa
n₀ = 71.6 # Limit adsoption / mol/kg

adsorption_params = HydrogenAdsorption.MDAParameters(n₀, p₀, α, β, m)

# Dubinin-Astakov parameters
P_lim = 77.75e6     # Pa
ψ = 7.3235          # mmol g-1
β = -0.0088         # mol kg-1 K-1
κ = 772.92          # J mol-1
γ = 18.828
m = 2.0            # Exponent in the isotherm equation
DA_params = HydrogenAdsorption.DAParameters(P_lim, ψ, β, κ, γ, m)