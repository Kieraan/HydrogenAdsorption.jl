using Revise
using HydrogenAdsorption
using Sundials
using Statistics
using Plots

# Parameters
# Material properties
# Activated carbon properties
cₛ = 825            # Specific heat capacity of carbon / J/kg K
mₛ = 0.440          # Mass of activated carbon / kg
ρₛ = 1990           # Density of activated carbon / kg/m³
kₛ = 0.9            # Thermal conductivity of activated carbon / W/m K [Not Necessary for this case]
k_eff = 0.21        # Effective thermal conductivity / W/(m·K)

# Hydrogen properties
cₚ = 14700.0        # Specific heat capacity of hydrogen / J/kg K
M_H2 = 2.0159e-3    # Molar mass of hydrogen / kg/mol
R = 8.314           # Ideal gas constant / J/mol K
k_g = 0.206         # Thermal conductivity of hydrogen / W/m K

# Generate struct for material properties
material_props = MaterialProperties(ρₛ, cₛ, mₛ, kₛ, ε_b, cₚ, M_H2, R, k_g, k_eff)

# Isotherm parameters
# Dubinin-Astakov parameters
P_lim = 77.75e6     # Pa
ψ = 7.3235          # mmol g-1
β = -0.0088         # mol kg-1 K-1
κ = 772.92          # J mol-1
γ = 18.828
m = 2.0            # Exponent in the isotherm equation
DA_params = DAParameters(P_lim, ψ, β, κ, γ, m)

# Tank parameters
L = 255e-3          # Length of the tank / m
R_T = 96e-3         # Interal radius of the tank / m
V = π * R_T^2 * L   # Volume of tank / m³
U = 36              # Heat transfer coefficient / W/(m²·K)
T₀ = 295.0          # Initial temperature of the tank / K
T_air = 295.0       # Ambient temperature / K
dr = 0.00025 / 2    # Radial step size / m
n_r, r_span, A, b = coefficient_matrix(R_T, dr, k_eff, U, T_air)
geometric_params = GeometricParameters(n_r, dr, V, A, b, r_span, R_T)

# Operational parameters
m_in = 9.4e-4 # Mass flow rate of hydrogen / m3/s
operational_params = OperationalParameters(U, T_air, m_in)

# Find initial conditions
Tᵢ = ones(n_r) * T₀ # Initial temperature / K
Pᵢ = 0.102564103e6 # Initial pressure / Pa
nₐᵢ = adsorption_isotherm(DA_params, Pᵢ, Tᵢ)
ρᵢ = ideal_gas_equation(Tᵢ, R, M_H2, P=Pᵢ)

# Find intiial conditions with new function
u₀, du₀, differential_vars = dae_setup(DA_params, material_props, geometric_params, operational_params, Pᵢ, T₀)
@assert u₀[1:n_r] == Tᵢ
@assert u₀[n_r+1:2*n_r] == nₐᵢ
@assert u₀[2*n_r+1] == ρᵢ[1]
@assert u₀[2*n_r+2] == Pᵢ
@assert u₀[2*n_r+3:end] == ρᵢ

println("End of assertions")