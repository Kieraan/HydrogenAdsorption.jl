using Revise
using HydrogenAdsorption
using Sundials
using Statistics
using Plots
using CSV
using DataFrames
VERBOSE = false # Set to true to print more information during the simulation

# parameters
# Tank parameters
V = 2.4946e-3 # Volume of tank / m³
L = 0.4 # Length of the tank / m
e = 0.050 - 0.0469 # Wall thickness / m
c_wall = 468 # Specific heat capacity of the tank wall / J/kg K
m_wall = 3.714 # Mass of the tank wall / kg
k_wall = 13 # Thermal conductivity of the tank wall / W/m K

# Inputs needed for the coefficient_matrix function
R_T = sqrt(V / (pi * L)) # Radius of the tank / m
dr = (1 / 2)^2 * 0.0025 / 2 # Radial step size / m
U = 36 # Heat transfer coefficient / W/(m²·K)

# Tank surface areas
Ai = 2 * pi * R_T * L # Inner surface area of the tank / m²
A0 = 2 * pi * (R_T + e) * L # Outer surface area of the tank / m²


T₀ = 281.0 # Initial temperature of the tank / K
T_air = 281.0 # Ambient temperature / K
# T_air = 282.5 # Ambient temperature / K

#T_H2 = 281.6 # Temperature of the incoming hydrogen gas / K
T_H2 = 297.6 # Data from finite element model / K


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
cₚ = 14700.0 # Specific heat capacity of hydrogen / J/kg K
cᵥ = 10134.0 # Specific heat capacity at constant volume / J/kg K
M_H2 = 2.0159e-3 # Molar mass of hydrogen / kg/mol
R = 8.314 # Ideal gas constant / J/mol K
k_g = 0.206 # Thermal conductivity of hydrogen / W/m K

k_eff = (kₛ * (1 - ε_b) + k_g * ε_b) # Effective thermal conductivity / W/(m·K)
n_r, r_span, A, b = coefficient_matrix(R_T, dr, k_eff, U, T_air)
material_props = MaterialProperties(ρₛ, cₛ, mₛ, kₛ, ε_b, cₚ, cᵥ, M_H2, R, k_g, k_eff)

# Geometric parameters
geometric_params = GeometricParameters(n_r, dr, V, L, A, b, r_span, R_T, e)

# Operational Parameters
U = 1 * 36.0 # Heat transfer coefficient / W/(m²·K)
m_in = 1 * 2.023e-5 # Mass flow rate of hydrogen / kg / s
operational_params = OperationalParameters(U, T_air, m_in, T_H2)

# Adsorption system parameters
par = AdsorptionParameters(MDA_params, material_props, geometric_params, operational_params)

# Find initial conditions
Tᵢ = ones(n_r) * T₀ # Initial temperature / K
Pᵢ = 0.033e6 # Initial pressure / Pa
nₐᵢ = adsorption_isotherm(MDA_params, Pᵢ, Tᵢ)
ρᵢ = ideal_gas_equation(Tᵢ, R, M_H2, P=Pᵢ)

# Find intiial conditions with new function
u₀, du₀, differential_vars = dae_setup(MDA_params, material_props, geometric_params, operational_params, Pᵢ, T₀)

function adsorption!(out, du, u, p, t)
    # Unpack parameters
    # Structs
    MDA_params = p.isotherm
    material_props = p.material
    geometric_params = p.geometric
    operational_params = p.operational

    # Scalars from structs
    # Material properties
    ρₛ = material_props.ρₛ
    ε_b = material_props.ε_b
    M_H2 = material_props.M_H2
    R = material_props.R
    k_eff = material_props.k_eff

    # Geometric parameters
    n_r = geometric_params.n_r
    dr = geometric_params.dr
    V = geometric_params.V
    r_span = geometric_params.r_span
    R_T = geometric_params.R_T

    # Operational parameters
    U = operational_params.U
    m_in = operational_params.m_in

    # Unpack state variables
    T = u[1:n_r]
    nₐ = u[n_r+1:2*n_r]
    ρ_avg = u[2*n_r+1]
    P = u[2*n_r+2]
    ρ = u[2*n_r+3:end]

    # Isosteric heat of adsorption
    dH = isosteric_heat_of_adsorption(MDA_params, P, T)
    #dH = isosteric_heat_of_adsorption(P, du[2*n_r+2], T, du[1:n_r], R)

    # Heat equation
    out[1:n_r] .= du[1:n_r] .- heat_equation(material_props, geometric_params, operational_params, u, du, dH)

    # Neumann BC time derivative for the tank centre 
    out[1] = du[1] - (4 * du[2] - du[3]) / 3

    # Robin BC time derivative to apply method of lines
    out[n_r] = du[n_r] - (4 * du[n_r-1] - du[n_r-2]) / (3 + 2 * U * dr / k_eff)

    # Adsorption isotherm
    out[n_r+1:2*n_r] .= nₐ .- adsorption_isotherm(MDA_params, P, T)

    # Macroscopic mass balance
    # mean(n_a .* r_span) / R computes the average adsorption of H2 
    dna_avg = sum((du[n_r+2:2*n_r] .* r_span[2:end] + du[n_r+1:2*n_r-1] .* r_span[1:end-1]) / 2 * (2 / R_T^2) * dr)
    # dna_avg_simple = mean(du[n_r+1:2*n_r])
    # println("dna_avg: $dna_avg, dna_avg_simple: $dna_avg_simple")
    out[2*n_r+1] = du[2*n_r+1] - (m_in / (V * ε_b) - ρₛ * (1 - ε_b) * M_H2 / ε_b * dna_avg)
    #out[2*n_r+1] = du[2*n_r+1] - (m_in / (V * ε_b) - ρₛ  * M_H2 / ε_b * dna_avg)

    # Ideal gas equation
    out[2*n_r+2] = du[2*n_r+2] - ideal_gas_equation(T, du[1:n_r], R, M_H2, R_T, r_span, ρ_avg, du[2*n_r+1])

    # Density of the gas
    out[2*n_r+3:end] = P .- ideal_gas_equation(T, R, M_H2, ρ=ρ)
end

t₀ = 0.0 # Initial time // s
t_f = 1042 # Final time // s
tspan = (t₀, t_f) # Time span for the simulation

# Create the DAE problem
prob = DAEProblem(adsorption!, du₀, u₀, tspan, p=par, differential_vars=differential_vars);
prob = remake(prob, p=par);
# @btime sol = solve(prob, IDA(linear_solver=:LapackDense), progress=true) # 838.288 s for dr = 0.000025
sol = solve(prob, IDA(linear_solver=:LapackDense), progress=true)
# Extract the solution
t = sol.t
T = [sol.u[i][1:n_r] for i in 1:length(sol.u)]
nₐ = [sol.u[i][n_r+1:2*n_r] for i in 1:length(sol.u)]
ρ_avg = [sol.u[i][2*n_r+1] for i in 1:length(sol.u)]
P = [sol.u[i][2*n_r+2] for i in 1:length(sol.u)]
ρ = [sol.u[i][2*n_r+3:end] for i in 1:length(sol.u)]

### Plotting ###
r_span = range(0, stop=R_T, length=n_r) # Generates radial nodes // m
#generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 200, :tab20b) # Generates the profiles plot for time_step = 200s
#generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 100, :tab20b) # Generates the profiles plot for time_step = 100s
generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 50, :tab20b) # Generates the profiles plot for time_step = 50s

middle_index = div(n_r, 2) # Index of the middle node
println("Middle index: $middle_index, final index: $(n_r)")


T_center = [T[i][1] for i in eachindex(T)] # Temperature at the tank center
T_middle = [T[i][middle_index] for i in eachindex(T)] # Temperature at the middle of the tank
T_radius = [T[i][end] for i in eachindex(T)] # Temperature at the tank wall


# Reconstruct final free and adsorber hydrogen mass
#  in the tank 
T_avg = [mean(T[i]) for i in eachindex(T)]
n_avg = [mean(nₐ[i]) for i in eachindex(nₐ)]
ρ_avg_nodes = [mean(ρ[i]) for i in eachindex(ρ)]

m_H2_gas = ρ_avg .* V * ε_b # Mass of hydrogen in the gas phase / kg
m_H2_ads = n_avg .* mₛ * M_H2 # Mass of hydrogen in the adsorbed phase / kg
m_H2_total = m_H2_gas .+ m_H2_ads # Total mass of hydrogen in the tank / kg 

if VERBOSE
    println("Final mass of hydrogen in the gas phase: $(m_H2_gas[end] * 1000) g")
    println("Final mass of hydrogen in the adsorbed phase: $(m_H2_ads[end] * 1000) g")
    println("Final total mass of hydrogen in the tank: $(m_H2_total[end] * 1000) g")
    println("Final average temperature in the tank: $(T_avg[end]) K")
    println("Final pressure in the tank: $(P[end] / 1e6) MPa")
    println("Final average gas density in the tank: $(ρ_avg_nodes[end]) kg/m³")
    println("Final average adsorption in the tank: $(n_avg[end]) mol/kg_ads")
end

# Compare with m_in * t_f
println("Total mass of hydrogen injected: $(m_in * t_f * 1000) g")

# Experimental data
# Temperature profiles
df_exp_middle_temp = CSV.read("inputs/xiao_middle_temp_finite_element.csv", DataFrame)
df_exp_center_temp = CSV.read("inputs/xiao_center_temp_finite_element.csv", DataFrame)
df_exp_wall_temp = CSV.read("inputs/xiao_wall_temp_finite_element.csv", DataFrame)

# Pressure profile
df_exp_pressure = CSV.read("inputs/xiao_pressure_finite_element.csv", DataFrame)

# Mass profiles
df_exp_gas_mass = CSV.read("inputs/xiao_gas_mass_finite_element.csv", DataFrame)
df_exp_adsorbed_mass = CSV.read("inputs/xiao_adsorbed_mass_finite_element.csv", DataFrame)
df_exp_total_mass = CSV.read("inputs/xiao_total_mass_finite_element.csv", DataFrame)

plt = plot(size=(1600, 900), margin=Plots.cm)
plt = plot(plt, xlabel="Time / s", ylabel="Temperature / K", title="Temperature Profiles")

# Middle temperature
plot!(plt, df_exp_middle_temp.Time, df_exp_middle_temp.Temperature, label="Experimental Middle Temperature", color=:black, marker=:circle, linestyle=:dash)
plot!(plt, t, T_middle, label="Simulated Middle Temperature", color=:black, lw=2)

# Center temperature
plot!(plt, df_exp_center_temp.Time, df_exp_center_temp.Temperature, label="Experimental Center Temperature", color=:blue, marker=:circle, linestyle=:dash)
plot!(plt, t, T_center, label="Simulated Center Temperature", color=:blue, lw=2)

# Wall temperature
plot!(plt, df_exp_wall_temp.Time, df_exp_wall_temp.Temperature, label="Experimental Wall Temperature", color=:red, marker=:circle, linestyle=:dash)
plot!(plt, t, T_radius, label="Simulated Wall Temperature", color=:red, lw=2)
display(plt)

# Pressure profile
plt = plot(df_exp_pressure.Time, df_exp_pressure.Pressure, xlabel="Time / s", ylabel="Pressure / MPa", title="Pressure Profile", label="Experimental Pressure", color=:black, marker=:circle, size=(800, 600), margin=Plots.cm)
plot!(plt, t, P ./ 1e6, label="Simulated Pressure", color=:black, lw=2)
display(plt)

# Mass profiles
plt = plot(size=(1600, 900), margin=Plots.cm)
plt = plot(plt, xlabel="Time / s", ylabel="Mass / g", title="Mass Profiles")
# Gas mass
plot!(plt, df_exp_gas_mass.Time, df_exp_gas_mass.Mass .* 1000, label="Experimental Gas Mass", color=:red, marker=:circle)
plot!(plt, t, m_H2_gas .* 1000, label="Simulated Gas Mass", color=:red, lw=2)
# Adsorbed mass
plot!(plt, df_exp_adsorbed_mass.Time, df_exp_adsorbed_mass.Mass .* 1000, label="Experimental Adsorbed Mass", color=:black, marker=:circle)
plot!(plt, t, m_H2_ads .* 1000, label="Simulated Adsorbed Mass", color=:black, lw=2)
# Total mass
plot!(plt, df_exp_total_mass.Time, df_exp_total_mass.Mass .* 1000, label="Experimental Total Mass", color=:blue, marker=:circle)
plot!(plt, t, m_H2_total .* 1000, label="Simulated Total Mass", color=:blue, lw=2)
display(plt)


