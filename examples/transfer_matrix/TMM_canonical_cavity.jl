using TransferMatrix
using GLMakie

λ_0 = 5.0
λs = range(3.2, 11, length = 700)
νs = 10^4 ./ λs
nperiods = 4
n1, n2 = 1.2, 2.6
t1 = λ_0 / (4 * n1)
t2 = λ_0 / (4 * n2)
t_middle = λ_0 / 2

air = Layer(λs, ones(length(λs)), zeros(length(λs)), t_middle)
l1 = Layer(λs, fill(n1, length(λs)), zeros(length(λs)), t1)
l2 = Layer(λs, fill(n2, length(λs)), zeros(length(λs)), t2)
dbr_unit = [l1, l2]
layers = [air, repeat(dbr_unit, nperiods)..., air, repeat(reverse(dbr_unit), nperiods)..., air]

# Tpp, Tss, Rpp, Rss = map(x -> collect(x), zip((calculate_tr(λ, layers) for λ in λs)...))


fig = Figure(size = (900, 600))
ax1 = Axis(fig[1, 1])
lines!(νs, Tpp, label = "pump off")

fig