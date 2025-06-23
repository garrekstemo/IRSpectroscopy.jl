using Revise
using MolecularPolaritonSpectra
using GLMakie

kB = 1.380649e-23  # Boltzmann constant in J/K
ħ = 1.0545718e-34  # Reduced Planck constant in J·s
c = 2.99792458e10  # Speed of light in cm/s

νs = 1900:0.1:2100

ν_c = ν_v = 2000.0  # in 1/cm
ω_v = ν_v * c  # convert to angular frequency in rad/s
ω_c = ν_c * c  # convert to angular frequency in rad/s
λ_v = 10^4 / ω_v

Δ = 20.0  # anharmonicity: ω_12 = ω_v - 2Δ
ν_12 = ν_v - 2 * Δ  # v=1 to v=2 transition

γ1 = 10.0  # fwhm of molecular vibrations induced by interaction with local bath
γ2 = 20.0
κ = 15.0
κL = κR = 0.5 * κ
N = 5e6

μ = 1e8
temp = 300.0  # temperature in Kelvin
β = 1 / (kB * temp)  # inverse temperature
g = abs(λ_v * μ) * tanh(β * ħ * ω_v / 2)

T_off = linear_transmission.(νs, ν_c, ν_v, N, g, κL, κR, κ, γ1)
T_on = linear_transmission_anharmonic.(νs, ν_c, ν_v, ν_12, N, g, 0.01, κL, κR, κ, γ1, γ2)
A = linear_absorption.(νs, ν_c, ν_v, κL, κ, γ1, N, g)
R = 1 .- T_off .- A
ΔT = T_on .- T_off
log_ΔT = -log.(T_on ./ T_off)


fig = Figure(size = (900, 800))
DataInspector()

ax1 = Axis(fig[1, 1], ylabel = "Transmission")
lines!(νs, T_off, label = "off")
lines!(νs, T_on, label = "on")
# lines!(νs, A, label = "Linear absorption")
# lines!(νs, R, label = "Linear reflection")

ax2 = Axis(fig[2, 1], xlabel = "Frequency (cm⁻¹)", ylabel = "ΔA")
# lines!(νs, -ΔT, label = "Linear response")
lines!(νs, log_ΔT, label = "Log-linear response")
axislegend(ax1)

fig