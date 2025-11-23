using Revise
using HydrogenAdsorption
using Sundials
using Statistics
using Plots
using CSV
using DataFrames

# --- Load Experimental Data for Full Cycle Plotting ---
pressure_exp_charge_data = CSV.read("inputs/xiao_pressure_finite_element.csv", DataFrame)
pressure_exp_dormant_data = CSV.read("inputs/xiao_pressure_finite_element_dormant.csv", DataFrame)
pressure_exp_discharge_data = CSV.read("inputs/xiao_pressure_finite_element_discharge.csv", DataFrame)

# --- Combine Experimental Data ---
pressure_full_cycle_exp = vcat(pressure_exp_charge_data, pressure_exp_dormant_data, pressure_exp_discharge_data)

# --- Load Time Series Data for Full Cycle Plotting ---
time_series_charge = CSV.read("outputs/charge_timeseries_data.csv", DataFrame)
time_series_dormant = CSV.read("outputs/dormant_timeseries_data.csv", DataFrame)
time_series_discharge = CSV.read("outputs/discharge_timeseries_data.csv", DataFrame)

# --- Combine Time Series Data ---
time_series_full_cycle = vcat(time_series_charge, time_series_dormant, time_series_discharge)
pressure_full_cycle = time_series_full_cycle[!, :Pressure_Pa]

# -----------------------------------------------------
# --- Load Temperature Data for Full Cycle Plotting ---
# -----------------------------------------------------

temperature_charge_data = CSV.read("outputs/charge_temperature_profiles.csv", DataFrame)
temperature_dormant_data = CSV.read("outputs/dormant_temperature_profiles.csv", DataFrame)
temperature_discharge_data = CSV.read("outputs/discharge_temperature_profiles.csv", DataFrame)

# --- Combine Temperature Data ---
temperature_full_cycle = hcat(temperature_charge_data, temperature_dormant_data[:, 3:end], temperature_discharge_data[:, 3:end])

# --- Extract Centre Temperature Data ---


# --- Extract Middle Node Temperature Data ---
middle_index = div(length(temperature_full_cycle[!, 1]), 2) # Index of the middle node


# --- Extract Wall Temperature Data ---


# --- 1. Global Plot Styling ---
# Apply a theme with larger fonts, similar to the Python 'seaborn-talk' style
theme(:default,
    fontfamily="sans-serif",
    titlefontsize=22,
    guidefontsize=20,     # X/Y axis labels
    tickfontsize=16,
    legendfontsize=10,
    margin=1.5Plots.cm,   # Add generous margins
    grid=true,
    gridstyle=:dash,
    gridalpha=0.5,
    framestyle=:box,      # Adds a full box, looks professional
    palette=:seaborn_pastel   # Use a modern, colorblind-friendly palette
)

# Define common line/marker styles for clarity
# (Simulated = thick solid, Experimental = thinner dash + markers)
sim_lw = 2.5
exp_lw = 1.5
exp_marker = (:circle, 6, 0.7) # (shape, size, alpha)


# Plot B: Pressure
p_pressure = plot(
    title="Pressure Profile",
    ylabel="Pressure / MPa",
    xlabel="Time / s",
    legend=:best, # Move legend inside
    legend_columns=2,     # Arrange in 2 columns (1 row)
    size=(1200, 900)
)

plot!(p_pressure, time_series_full_cycle[!, :Time_s], pressure_full_cycle ./ 1e6,
    label="Simulated", color=1, lw=sim_lw)
plot!(p_pressure, pressure_full_cycle_exp.Time, pressure_full_cycle_exp.Pressure,
    label="Experimental", color=1, marker=exp_marker, linestyle=:dash, lw=exp_lw)

display(p_pressure)