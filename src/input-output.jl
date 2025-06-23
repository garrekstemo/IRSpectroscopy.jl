"""
    b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, A)
"""
function b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, f_pu, Δ, γ_m, g_1, g_3)

    A_ω = A(ω, f_pu, ω_0, Δ, γ_m, g_1, g_3)
    -im * 0.5 * κ * (ω - ω_0 + im * γ_1) * (ω - ω_12 + im * γ_3) / ( (ω - ω_c + im * 0.5 * κ) * (ω - ω_0 + im * γ_1) * (ω - ω_12 + im * γ_3) - N * A_ω)
end

"""
    A(ω, f_pu, ω_0, Δ, γ_m, g_1, B)
"""
function A(ω, f_pu, ω_0, Δ, γ_m, g_1, g_3)
    g_1^2 * (ω - (ω_0 - 2*Δ) + 3 * 0.5 * im * γ_m) + f_pu * B(ω, ω_0, Δ, γ_m, g_1, g_3)
end

"""
    B(ω, ω_0, Δ, γ_m, g_1, g_3)

g_3 = c * g_1, where c = g_3/g_1 can be positive or negative.
Here we use c.
"""
function B(ω, ω_0, Δ, γ_m, g_1, g_3)
    g_3 * (2 * (2 * g_1 + g_3) * (ω - ω_0) + (4 * g_1 + g_3) * im * γ_m) - 4 * Δ * g_1^2
end

"""
    io_transmission(ω, f_pu, ω_0, ω_c, Δ, γ_m, κ, g_1, g_3, N)

Single-pulse transmission spectrum.

ω = incident wavelength
f_pu = fraction of molecules excited by pump pulse (0 = no pumping)
ω_0 = molecular fundamental frequency
ω_c = cavity mode frequency
Δ = molecular anharmonicity, Δ = ω_0 - ω_12 (always positive)
γ_m = molecular mode line width (fwhm)
κ/2 = cavity line width (but papers seem to maybe just use κ?)
g_1 = single-molecule light-matter coupling
g_3 = deviation of vibration 1->2 transition dipole moment μ_12 
      from the corresponding harmonic oscillator:

      μ_12 = √2 * μ_1 * (1 + g_3 / g_1)
"""
function io_transmission(ω, f_pu, ω_0, ω_c, Δ, γ_m, κ, g_1, g_3, N)

    ω_12 = ω_0 - 2 * Δ
    γ_1 = 0.5 * γ_m
    γ_3 = 0.5 * 3 * γ_m

    A_ω = A.(ω, f_pu, ω_0, Δ, γ_m, g_1, g_3)
    T = abs2.(b_out.(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, f_pu, Δ, γ_m, g_1, g_3))

    return T
end
