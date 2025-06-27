using MolecularPolaritonSpectra
using GLMakie

ω = 1900:0.1:2100

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
Δ = 22.0  # anharmonicity: ω_12 = ω_v - 2Δ
ω_12 = ω_v - 2 * Δ  # v=1 to v=2 transition

T_off = io_transmission(ω, 0, ω_v, ω_c, Δ, γ_m, κ, g_1, g_3, N)
T_on = io_transmission(ω, f_pu, ω_v, ω_c, Δ, γ_m, κ, g_1, g_3, N)
ΔT = T_on .- T_off
# ΔT = T_on ./ T_off
ΔA = -log10.(T_on ./ T_off)


fig = Figure(size = (900, 800))

ax1 = Axis(fig[1, 1], ylabel = "Transmission")
lines!(ω, T_off, label = "pump off")
lines!(ω, T_on, label = "pump on")

ax2 = Axis(fig[2, 1], xlabel = "Frequency (cm⁻¹)", ylabel = "ΔA")
# lines!(ω, -ΔT, label = "Nonlinear response")
lines!(ω, ΔA, label = "Log-linear response")
# vlines!(ω_12, linestyle = :dash, color = :firebrick4, label = "Molecular excited state")

fig