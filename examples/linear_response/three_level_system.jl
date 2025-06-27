using Revise
using MolecularPolaritonSpectra
using GLMakie

kB = 1.380649e-23  # Boltzmann constant in J/K
ħ = 1.0545718e-34  # Reduced Planck constant in J·s
c = 2.99792458e10  # Speed of light in cm/s

ωs = -5:0.01:5

ω_12 = ω_ph = 1
ω_23 = 2ω_ph
κ = 0.1
γ = 0.3

μ = 1e8
temp = 300.0  # temperature in Kelvin
β = 1 / (kB * temp)  # inverse temperature
N = 1
g = 1 #abs(λ_v * μ) * tanh(β * ħ * ω_v / 2)

κL = κR = 1 # 0.5 * κ

susceptibility = χ.(ωs, ω_12, N, g, γ) .+ χ.(ωs, ω_23, N, g, γ)
T = linear_transmission(ωs, susceptibility, ω_ph, κL, κR, κ)
# A = linear_absorption.(ωs, ω_ph, ω_12, κL, κ, γ1, N, g)
# R = 1 .- T .- A


fig = Figure(size = (900, 800))
DataInspector()

ax1 = Axis(fig[1, 1], ylabel = "Intensity")
lines!(ωs, T, label = "Transmission")

ax2 = Axis(fig[2, 1], xlabel = "Frequency (cm⁻¹)", ylabel = "Intensity")

fig