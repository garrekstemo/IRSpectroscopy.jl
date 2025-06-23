using Revise
using IRSpectroscopy
using GLMakie

ω = 1900:0.1:2100

# ---------------------------------------------------------------------- #
# Input-Output Theory

γ_m = 10.0  # fwhm of molecular vibrations induced by interaction with local bath
γ_1 = 0.5 * γ_m
γ_3 = 3 * γ_m / 2
κ = 10.0
N = 1
f_pu = 0.01

ω_v = 2000.0
ω_c = 2000.0
g_1 = 25 / √(N)
g_3 = 0.2
Δ = 20.0  # anharmonicity: ω_12 = ω_v - 2Δ
ω_12 = ω_v - 2 * Δ  # v=1 to v=2 transition

T_off = io_transmission(ω, 0, ω_v, ω_c, Δ, γ_m, κ, g_1, g_3, N)
T_on = io_transmission(ω, f_pu, ω_v, ω_c, Δ, γ_m, κ, g_1, g_3, N)
ΔT = T_on .- T_off

# ---------------------------------------------------------------------- #
# Linear Response Theory

Tlin_off = linear_transmission.(ω, ω_c, ω_v, N, g_1, 0.5 * κ, 0.5 * κ, κ, γ_m)
Tlin_on = linear_transmission_anharmonic.(ω, ω_c, ω_v, ω_12, N, g_1, 0.02, 0.5 * κ, 0.5 * κ, κ, γ_m, γ_3)
Alin = linear_absorption.(ω, ω_c, ω_v, 0.5 * κ, κ, γ_m, N, g_1)
Rlin = 1 .- Tlin_off .- Alin
ΔTlin = Tlin_on .- Tlin_off

# ---------------------------------------------------------------------- #


fig = Figure(size = (900, 800))

ax1 = Axis(fig[1, 1])
lines!(ω, T_off, label = "pump off")
lines!(ω, T_on, label = "pump on")

ax2 = Axis(fig[2, 1])
lines!(ω, ΔT, label = "Nonlinear response")
# vlines!(ω_12, linestyle = :dash, color = :firebrick4, label = "Molecular excited state")

ax3 = Axis(fig[1, 2])
lines!(ω, Tlin_off, label = "Linear response")
lines!(ω, Tlin_on, label = "Linear response")
# lines!(ω, Alin, label = "Linear absorption")
# lines!(ω, Rlin, label = "Linear reflection")

ax4 = Axis(fig[2, 2])
lines!(ω, ΔTlin, label = "Linear response")
axislegend(ax1)



fig