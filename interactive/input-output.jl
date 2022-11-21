using Revise
using GLMakie
using IRSpectroscopy
using SomeMakieThemes

set_theme!(theme_retina())

ω = 2050:0.1:2300

γ_m = 10.0  # fwhm of molecular vibrations induced by interaction with local bath
γ_1 = 0.5 * γ_m
γ_3 = 3 * γ_m / 2
κ = 10.0
N = 1
f_pu = 0.1

ω_0 = Observable(2170.0)   # fundamental vibrational frequency
ω_c = Observable(2170.0)  # cavity mode frequency
ω_12 = Observable(2150.0)  # v=1 to v=2 transition
g_1 = Observable(25 / √(N))
Δ = Observable((to_value(ω_0) - to_value(ω_12)) / 2)  # anharmonicity: ω_12 = ω_0 - 2Δ
g_3 = Observable(0.2)

T0 = transmission(ω, 0, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
T_fpu = transmission(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
ΔT = pp_spectrum(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)

oToff = Observable(T0)
oTon = Observable(T_fpu)
oΔT = Observable(ΔT)


fig = Figure(figure_padding = 20, resolution = (2300, 1900))
display(fig)

Label(fig[0, :], "Vibrational Polariton Nonlinear Spectroscopy")

tbΔ = Textbox(fig, placeholder= "Set anharmonicity...", width = 400)
tbω_0 = Textbox(fig, placeholder = "$(round(Int, to_value(ω_0)))cm⁻¹", width = 300)

fig[2, 1][1, 1] = vgrid!(
    Label(fig, "ω_0 (cm⁻¹)"),
    Label(fig, "Δ (cm⁻¹)");
    tellwidth = false,
    tellheight = false,
)
fig[2, 1][1, 2] = vgrid!(
    tbω_0,
    tbΔ;
    tellwidth = false,
    tellheight = false,
)

sg = SliderGrid(
    fig[2, 1][2, 1:2],
    (label = L"ω_c", range = 2100:1:2300, format = "{:} cm⁻¹", startvalue = 2170),
    (label = L"g_1", range = 5:0.1:50, format = "{:.1f} cm⁻¹", startvalue = 20.0),
    (label = L"\frac{g_3}{g_1}", range = -5:0.1:5, format = "{:.1f}", startvalue = 0.0),
    (label = L"Δ", range = 0:0.1:20, format = "{:.1f} cm⁻¹", startvalue = 0.0),
    tellheight = false,
    tellwidth = false
)

ax = Axis(fig[1, 1:2], title = "Pump-probe spectrum", xlabel = "ω (cm⁻¹)", ylabel = L"ΔT = T_{f^{pu}} - T_0 (\text{a.u.})", 
            xticks = LinearTicks(10), yticks = LinearTicks(5),
            yticklabelspace = 100.0  # plot won't resize when tick digits increase/decrease
            )

lines!(ax, ω, oΔT)
vlines!(ax, ω_12, linestyle = :dash, color = :firebrick4)

ax2 = Axis(fig[2, 2], title = "Transmission spectrum", xlabel = "ω (cm⁻¹)", ylabel = "T (arb.)", 
            xticks = LinearTicks(5), yticks = LinearTicks(5))

lines!(ax2, ω, oToff, label = "pump off")
lines!(ax2, ω, oTon, label = "pump on")
axislegend(ax2)

sliderobservables = [s.value for s in sg.sliders]
connect!(ω_c, sliderobservables[1])
connect!(g_1, sliderobservables[2])
connect!(g_3, sliderobservables[3])
connect!(Δ, sliderobservables[4])

on(tbω_0.stored_string) do s
    ω_0[] = parse(Float64, s)
    oToff[] = transmission(ω, 0, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oTon[] = transmission(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oΔT[] = pp_spectrum(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    autolimits!(ax)
    autolimits!(ax2)
end
on(tbΔ.stored_string) do s
    Δ[] = parse(Float64, s)
    oToff[] = transmission(ω, 0, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oTon[] = transmission(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oΔT[] = pp_spectrum(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    autolimits!(ax)
    autolimits!(ax2)
end
on(ω_c) do val
    # Δ[] = abs(to_value(ω_c) - to_value(ω_0)) / 2
    oToff[] = transmission(ω, 0, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oTon[] = transmission(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oΔT[] = pp_spectrum(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    autolimits!(ax)
    autolimits!(ax2)
end 

on(g_1) do val
    oToff[] = transmission(ω, 0, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oTon[] = transmission(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oΔT[] = pp_spectrum(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    autolimits!(ax)
    autolimits!(ax2)
end

on(g_3) do val
    oToff[] = transmission(ω, 0, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oTon[] = transmission(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oΔT[] = pp_spectrum(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    autolimits!(ax)
    autolimits!(ax2)
end

on(Δ) do val
    ω_12[] = to_value(ω_0) - 2 * to_value(Δ)
    oToff[] = transmission(ω, 0, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oTon[] = transmission(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    oΔT[] = pp_spectrum(ω, f_pu, to_value(ω_0), to_value(ω_c), to_value(Δ), γ_m, κ, to_value(g_1), to_value(g_3), N)
    autolimits!(ax)
    autolimits!(ax2)
end
