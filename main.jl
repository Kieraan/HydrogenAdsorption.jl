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