using GLMakie
using IRSpectroscopy

ω = 2050:0.1:2300

γ_m = 10.0  # fwhm of molecular vibrations induced by interaction with local bath
γ_1 = 0.5 * γ_m
γ_3 = 3 * γ_m / 2
κ = 10.0
N = 1
f_pu = 0.1

ω_0 = Observable(2170.0)   # fundamental vibrational frequency
ω_c = Observable(2170.0)  # cavity mode frequency
g_1 = Observable(25 / √(N))
g_3 = Observable(0.2)
Δ = Observable(20.0)  # anharmonicity: ω_12 = ω_0 - 2Δ
ω_12 = @lift($ω_0 - 2 * $Δ)  # v=1 to v=2 transition

Toff = @lift(transmission(ω, 0, $(ω_0), $(ω_c), $(Δ), γ_m, κ, $(g_1), $(g_3), N))
Ton = @lift(transmission(ω, f_pu, $(ω_0), $(ω_c), $(Δ), γ_m, κ, $(g_1), $(g_3), N))
ΔT = @lift(pp_spectrum(ω, f_pu, $(ω_0), $(ω_c), $(Δ), γ_m, κ, $(g_1), $(g_3), N))


fig = Figure(size = (850, 800))

tbω_0 = Textbox(fig, placeholder = "$(round(Int, to_value(ω_0))) cm⁻¹", width = 300)
flipbutton = Button(fig, label = "Flip ΔT", width = 200)

fig[2, 1][1, 1] = vgrid!(
    Label(fig, "Fundamental frequency"),    
    tbω_0,
    flipbutton,
    tellwidth = false,
    tellheight = false,
)

sg = SliderGrid(
    fig[2, 1][2, 1],
    (label = L"ω_c", range = 2100:1:2300, format = "{:} cm⁻¹", startvalue = to_value(ω_0)),
    (label = L"g_1", range = 5:0.1:50, format = "{:.1f} cm⁻¹", startvalue = 20.0),
    (label = L"\frac{g_3}{g_1}", range = -5:0.1:5, format = "{:.1f}", startvalue = 0.0),
    (label = L"Δ", range = 0:0.1:20, format = "{:.1f} cm⁻¹", startvalue = 10.0),
    tellheight = false,
    tellwidth = false
)

ax1 = Axis(fig[1, 1:2], title = "Pump-probe spectrum", xlabel = "ω (cm⁻¹)", ylabel = L"ΔT = T_{f^{pu}} - T_0 \text{ (arb.)}", 
            xticks = LinearTicks(10), yticks = LinearTicks(5),
            yticklabelspace = 50.0  # plot won't resize when tick digits increase/decrease
            )

lines!(ax1, ω, ΔT, label = "Nonlinear response")
vlines!(ax1, ω_12, linestyle = :dash, color = :firebrick4, label = "Molecular excited state")
axislegend(ax1)

ax2 = Axis(fig[2, 2], title = "Transmission spectrum", xlabel = "ω (cm⁻¹)", ylabel = L"T \text{ (arb.)}", 
            xticks = LinearTicks(5), yticks = LinearTicks(5))

lines!(ax2, ω, Toff, label = "pump off")
lines!(ax2, ω, Ton, label = "pump on")
axislegend(ax2)

Label(fig[0, :], "\nQuantum model for cavity pump-probe spectroscopy\ninput-output theory\n ")


sliderobservables = [s.value for s in sg.sliders]
connect!(ω_c, sliderobservables[1])
connect!(g_1, sliderobservables[2])
connect!(g_3, sliderobservables[3])
connect!(Δ, sliderobservables[4])

for slider in sg.sliders
    on(slider.value) do _
        autolimits!(ax1)
        autolimits!(ax2)
    end
end

on(tbω_0.stored_string) do s
    ω_0[] = parse(Float64, s)
end

on(flipbutton.clicks) do _
    ΔT[] *= -1
    autolimits!(ax1)
end

fig