using Revise
using TransferMatrix
using MolecularPolaritonSpectra

λ_0 = 5.0
λs = range(4.76, 5.26, step=0.002)
νs = 10^4 ./ λs
nperiods = 5
n1, n2 = 1.3, 2.2
t1 = λ_0 / (4 * n1)
t2 = λ_0 / (4 * n2)

A = 5e3
excited_fraction = 0.015
A1 = A * (1 - excited_fraction)
A2 = A * excited_fraction
ν_0 = 10^4 / λ_0
ν_12 = 1980
Γ = 10.0
Γ_12  = 20.0

n_bg = 1.0
n, k = nk(νs, A, ν_0, Γ, n_bg)
n12, k12 = nk2(νs, (A1, ν_0, Γ, A2, ν_12, Γ_12, n_bg))
t_middle = λ_0 / (2 * n_bg) .- 0.0029

air = Layer(λs, ones(length(λs)), zeros(length(λs)), t_middle)
l1 = Layer(λs, fill(n1, length(λs)), zeros(length(λs)), t1)
l2 = Layer(λs, fill(n2, length(λs)), zeros(length(λs)), t2)
l_middle = Layer(λs, n, k, t_middle)
l_middle2 = Layer(λs, n12, k12, t_middle)
dbr_unit = [l1, l2]
layers = [air, repeat(dbr_unit, nperiods)..., l_middle, repeat(reverse(dbr_unit), nperiods)..., air]
layers2 = [air, repeat(dbr_unit, nperiods)..., l_middle2, repeat(reverse(dbr_unit), nperiods)..., air]

Tpp_off, Tss, Rpp, Rss = map(x -> collect(x), zip((calculate_tr(λ, layers) for λ in λs)...))
Tpp_on, Tss_on, Rpp_on, Rss_on = map(x -> collect(x), zip((calculate_tr(λ, layers2) for λ in λs)...))


fig = Figure(size = (800, 600))
DataInspector()
ax1 = Axis(fig[1, 1], xlabel = "Wavelength (nm)", ylabel = "Transmission")
lines!(ax1, λs, Tpp_off, label = "Tpp off")
lines!(ax1, λs, Tss, label = "Tss off")

ax2 = Axis(fig[2, 1], xlabel = "Wavelength (nm)", ylabel = "Reflection")
lines!(νs, Tpp_off .- Tpp_on)
fig