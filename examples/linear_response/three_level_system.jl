using Revise
using MolecularPolaritonSpectra
using GLMakie

ωs = -5:0.01:7
ω_ph = ω_12 = 1
ω_23 = 2ω_ph
κ = 0.1
κL = κR = 0.5 * κ
γ = 0.3

N = 1
g = 1

μ_zy = [
    0 0 0
    1 0 0
    1 1 0
]
ω_zy = [
    0    0    0
    ω_12 0    0
    3    ω_23 0
]
ps = [
    [1, 1, 1],
    [0.45, 0.45, 0.1],
    [0.7, 0.2, 0.1],
]

fig = Figure(size = (600, 1100))
DataInspector()

axes = [Axis(fig[i, 1], ylabel = "Intensity") for i in eachindex(ps)]
for i in eachindex(ps)
    susceptibility = χ(ωs, ps[i], μ_zy, ω_zy, N, γ, 1)
    T = linear_transmission(ωs, susceptibility, ω_ph, κL, κR, κ)
    A = ω_ph / κR * imag.(susceptibility) .* T
    R = 1 .- T .- A
    lines!(axes[i], ωs, T, label = "T")
    lines!(axes[i], ωs, A, label = "A")
    lines!(axes[i], ωs, R, label = "R")
    
    axislegend(axes[i], position = :rc)
end

fig
