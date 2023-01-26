"""
    b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, A)
"""
function b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, A)

    bout = Vector{ComplexF64}(undef, length(ω))

    for (i, ω_i) in enumerate(ω)
        bout_i = -im * 0.5 * κ * (ω_i - ω_0 + im * γ_1) * (ω_i - ω_12 + im * γ_3) / ( (ω_i - ω_c + im * 0.5 * κ) * (ω_i - ω_0 + im * γ_1) * (ω_i - ω_12 + im * γ_3) - N * A[i] )
        bout[i] = bout_i
    end

    return bout
end

"""
    A(ω, f_pu, ω_0, Δ, γ_m, g_1, B)
"""
function A(ω, f_pu, ω_0, Δ, γ_m, g_1, B)

    a = Vector{ComplexF64}(undef, length(ω))

    for (i, ω_i) in enumerate(ω)
        a_i = g_1^2 * (ω_i - (ω_0 - 2*Δ) + 3 * 0.5 * im * γ_m) + f_pu * B[i]
        a[i] = a_i
    end
    return a
end

"""
    B(ω, ω_0, Δ, γ_m, g_1, g_3_ratio)

g_3 = c * g_1, where c = g_3/g_1 can be positive or negative.
Here we use c.
"""
function B(ω, ω_0, Δ, γ_m, g_1, g_3_ratio)

    b = Vector{ComplexF64}(undef, length(ω))

    g_3 = g_3_ratio * g_1
    for (i, ω_i) in enumerate(ω)
        b_i = g_3 * (2 * (2 * g_1 + g_3) * (ω_i - ω_0) + (4 * g_1 + g_3) * im * γ_m) - 4 * Δ * g_1^2
        b[i] = b_i
    end
    b
end

"""
    transmission(ω, f_pu, ω_0, ω_c, Δ, γ_m, κ, g_1, g_3_ratio, N)

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
      
      (we use just the ratio, g_3 / g_1)
"""
function transmission(ω, f_pu, ω_0, ω_c, Δ, γ_m, κ, g_1, g_3_ratio, N)

    ω_12 = ω_0 - 2 * Δ
    γ_1 = 0.5 * γ_m
    γ_3 = 0.5 * 3 * γ_m

    B_ω = B(ω, ω_0, Δ, γ_m, g_1, g_3_ratio)
    A_ω = A(ω, f_pu, ω_0, Δ, γ_m, g_1, B_ω)
    T = abs2.(b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, A_ω))

    return T
end

"""
    pp_spectrum(ω, f_pu, ω_0, ω_c, Δ, γ_1, κ, g_1, g_3_ratio, N)

Pump-probe spectrum
    ΔT = T_fpu - T_0
"""
function pp_spectrum(ω, f_pu, ω_0, ω_c, Δ, γ_m, κ, g_1, g_3_ratio, N)

    ω_12 = ω_0 - 2 * Δ
    γ_1 = 0.5 * γ_m
    γ_3 = 0.5 * 3 * γ_m

    B_ω = B(ω, ω_0, Δ, γ_m, g_1, g_3_ratio)
    A_ω = A(ω, f_pu, ω_0, Δ, γ_m, g_1, B_ω)
    T_fpu = abs2.(b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, A_ω))

    A0 = A(ω, 0.0, ω_0, Δ, γ_m, g_1, B_ω)
    T_0 = abs2.(b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, A0))

    return T_fpu .- T_0
end


"""
    pp_spectrum(ω, f_pu, ω_0, ω_c, Δ, γ_1, κ, g_1, g_3_ratio, N)

Pump-probe spectrum
    ΔT = T_fpu - T_0
"""
function pp_spectrum(ω, f_pu, ω_0, ω_c1, ω_c2, Δ, γ_m, κ, g_1, g_3_ratio, N)

    ω_12 = ω_0 - 2 * Δ
    γ_1 = 0.5 * γ_m
    γ_3 = 0.5 * 3 * γ_m

    B_ω = B(ω, ω_0, Δ, γ_m, g_1, g_3_ratio)
    A_ω = A(ω, f_pu, ω_0, Δ, γ_m, g_1, B_ω)
    T_fpu = abs2.(b_out(ω, ω_0, ω_c1, ω_12, γ_1, γ_3, κ, N, A_ω))

    A0 = A(ω, 0.0, ω_0, Δ, γ_m, g_1, B_ω)
    T_0 = abs2.(b_out(ω, ω_0, ω_c2, ω_12, γ_1, γ_3, κ, N, A0))

    return T_fpu .- T_0
end
