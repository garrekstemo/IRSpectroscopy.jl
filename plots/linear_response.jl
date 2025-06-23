using Revise
using MolecularPolaritonSpectra
using GLMakie

kB = 1.380649e-23  # Boltzmann constant in J/K
ħ = 1.0545718e-34  # Reduced Planck constant in J·s

ωs = 1900:0.1:2100

ω_c = ω_v = 2000.0
λ_v = 10^4 / ω_v 
Δ = 20.0  # anharmonicity: ω_12 = ω_v - 2Δ
ω_12 = ω_v - 2 * Δ  # v=1 to v=2 transition

γ1 = 10.0  # fwhm of molecular vibrations induced by interaction with local bath
γ2 = 20.0
κ = 10.0
κL = κR = 0.5 * κ
N = 1e3

μ = 1e10
temp = 300.0  # temperature in Kelvin
β = 1 / (kB * temp)  # inverse temperature
g = abs(λ_v * μ) * tanh(β * ħ * ω_v / 2)

T_off = linear_transmission.(ωs, ω_c, ω_v, N, g, κL, κR, κ, γ1)
T_on = linear_transmission_anharmonic.(ωs, ω_c, ω_v, ω_12, N, g, 0.01, κL, κR, κ, γ1, γ2)
A = linear_absorption.(ωs, ω_c, ω_v, κL, κ, γ1, N, g)
R = 1 .- T_off .- A
ΔT = T_on .- T_off


fig = Figure(size = (900, 800))

ax1 = Axis(fig[1, 1])
lines!(ωs, T_off, label = "off")
lines!(ωs, T_on, label = "on")
# lines!(ωs, A, label = "Linear absorption")
# lines!(ωs, R, label = "Linear reflection")

ax2 = Axis(fig[2, 1])
lines!(ωs, ΔT, label = "Linear response")
axislegend(ax1)

fig