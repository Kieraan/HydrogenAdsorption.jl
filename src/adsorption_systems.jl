"""
System of differential equations for the adsorption process in a hydrogen tank. 
This function computes the time derivative of the state variables, which include temperature, adsorption amount, average density, pressure, and local densities.
It uses the method of lines to discretize the spatial domain and solve the system of equations.

This system in particular uses the Modified Dubinin-Astakov isotherm for adsorption, which is suitable for hydrogen storage applications.

Inputs
- `out`: Output array
- `du`: Array containing the current time derivatives of the state variables.
- `u`: Current state of the system, including temperature, adsorption amount, average density, pressure, and local densities.
- `p`: Parameters for the system
- `t`: Current time (not used in this function, but included for compatibility with ODE solvers).

Outputs: The function modifies the `out` array in place

Description of variables:
- `T`: Array of dimension `n_r` representing the temperature at each radial node
- `nₐ`: Array of dimension `n_r` representing the adsorption amount at each radial node
- `ρ_avg`: Average density of hydrogen in the tank
- `P`: Pressure in the tank
- `ρ`: Array of dimension `n_r` representing the local densities of hydrogen at each radial node
"""
function MDA_adsorption!(out, du, u, p, t)
    # Unpacking variables
    T = u[1:n_r]
    nₐ = u[n_r+1:2*n_r]
    ρ_avg = u[2*n_r+1]
    P = u[2*n_r+2]
    ρ = u[2*n_r+3:end]

    dH = α .* ((log.(n₀ ./ nₐ)) .^ (1 / m))

    # Heat diffusion, implicit form
    # f(x) = 0 = dT/dt - A*T
    # A is the discretisation matrix
    # b should be zero as it is recalculated in 63.
    out[1:n_r] .= du[1:n_r] .- 1 / (ρₛ * cₛ * (1 - ε_b) + ρ_avg * cₚ * ε_b) * (A * T .+ dH .* du[n_r+1:2*n_r] ./ V .+ du[2*n_r+2])

    # Neumann BC time derivative for the tank centre 
    out[1] = du[1] - (4 * du[2] - du[3]) / 3

    # Robin BC time derivative to apply method of lines
    out[n_r] = du[n_r] - (4 * du[n_r-1] - du[n_r-2]) / (3 + 2 * U * dr / k_eff)

    # Adsorption equation
    out[n_r+1:2*n_r] .= nₐ .- n₀ .* exp.(-(R .* T ./ (α .+ β .* T)) .^ m .* (log.(p₀ / P)) .^ m)

    # Macroscopic mass balance 
    # mean(n_a .* r_span) / R computes the average adsorption of H2 
    out[2*n_r+1] = du[2*n_r+1] - (m_in / (V * ε_b) - ρₛ * (1 - ε_b) * M_H2 / ε_b * mean(du[n_r+1:2*n_r] .* r_span) / R_T)

    # Ideal gas pressure equation
    # mean(T .* r_span) / R_T computes the pressure inside the tank
    # out[2*n_r+2] = P - ρ_avg * R .* (mean(T .* r_span) / R_T) / M_H2
    # out[2*n_r+2] = du[2*n_r+2] - dPdt_fixed
    out[2*n_r+2] = du[2*n_r+2] - R / M_H2 * (du[2*n_r+1] * (mean(T .* r_span) / R_T) + ρ_avg * (mean(du[1:n_r] .* r_span) / R_T))

    # Pressure equation used to compute ρ in each node
    out[2*n_r+3:end] = P .- ρ .* R .* T ./ M_H2
end