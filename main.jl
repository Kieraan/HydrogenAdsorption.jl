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

n_r, r_span, A, Tᵢ, b = coefficient_matrix(R_T, dr, k_eff, U, T₀, T_air)

# Isotherm parameters
# Modified Dubinin-Astakov Isotherm parameters
α = 3080.0 # Enthalpic Factor / J/mol
β = 18.9 # Entropic Factor / J/mol K
m = 2.0 # Exponential factor
p₀ = 1470e6 # Saturation pressure / Pa
n₀ = 71.6 # Limit adsoption / mol/kg

MDA_params = MDAParameters(n₀, p₀, α, β, m)

# Dubinin-Astakov parameters
P_lim = 77.75e6     # Pa
ψ = 7.3235          # mmol g-1
β = -0.0088         # mol kg-1 K-1
κ = 772.92          # J mol-1
γ = 18.828
m = 2.0            # Exponent in the isotherm equation
DA_params = DAParameters(P_lim, ψ, β, κ, γ, m)

# Material properties
# Activated carbon properties
ρₛ = 517.6 # Density of activated carbon / kg/m³
cₛ = 825 # Specific heat capacity of carbon / J/kg K
mₛ = 0.671 # Mass of activated carbon / kg
kₛ = 0.646 # Thermal conductivity of activated carbon / W/m K
ε_b = 0.49 # Bed porosity

# Hydrogen properties
cₚ = 14700 # Specific heat capacity of hydrogen / J/kg K
M_H2 = 2.0159e-3 # Molar mass of hydrogen / kg/mol
R = 8.314 # Ideal gas constant / J/mol K
k_g = 0.206 # Thermal conductivity of hydrogen / W/m K

k_eff = (kₛ * (1 - ε_b) + k_g * ε_b) # Effective thermal conductivity / W/(m·K)
material_props = MaterialProperties(ρₛ, cₛ, mₛ, kₛ, ε_b, cₚ, M_H2, R, k_g, k_eff)

# Geometric parameters
geometric_params = GeometricParameters(n_r, dr, V, A, b, r_span, R_T)

# Operational Parameters
U = 36.0 # Heat transfer coefficient / W/(m²·K)
m_in = 2.023e-5 # Mass flow rate of hydrogen / kg / s
operational_params = OperationalParameters(U, m_in)

# Adsorption system parameters
p = AdsorptionParameters(MDA_params, material_props, geometric_params, operational_params)
println(typeof(p))