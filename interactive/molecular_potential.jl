using Revise
using GLMakie
using IRSpectroscopy
using SomeMakieThemes
set_theme!(theme_retina())

q = 0:0.001:40
q_eq = 10
De = 9.1
α = 0.2
ns = 0:6
ν0 = 2
ν_exe = 0.03

morse = morse_potential(q, q_eq, α, De)
harm = harmonic(q, q_eq, α, De)
levels = energy_levels(ns, ν0, ν_exe)
harm_levels = ν0 * (ns .+ 0.5)

f = Figure()
display(f)
DataInspector(f)

ax = Axis(f[1, 1], title = "Morse", xlabel = "Intermolecular distance")
lines!(q, morse, label = "Morse")
lines!(q, harm, label = "harmonic")
hlines!(levels, xmin=0.1, xmax=0.5)
hlines!(harm_levels, xmin = 0.1, xmax=0.5)
ylims!(-1, 20)
axislegend(ax)