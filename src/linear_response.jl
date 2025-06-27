function χ(ω, ω_v, N, g, γ)
    -N * abs2(g) / (ω - ω_v + im * 0.5 * γ)
end

"""
    linear_transmission(ω, ω_ph, ω_v, N, g, κL, κR, κ, γ)
Calculates the linear transmission spectrum for a system with a cavity mode and molecular vibrations.
ω: Frequency at which to evaluate the transmission.
ω_ph: Frequency of the cavity mode.
ω_v: Frequency of the molecular vibration.
N: Number of molecules.
g: Coupling strength between the cavity mode and the molecular vibration.
κL: Left cavity decay rate.
κR: Right cavity decay rate.
κ: Total cavity decay rate.
"""
function linear_transmission(ω, χ, ω_ph, κL, κR, κ)
    T = similar(ω)
    for i in eachindex(ω)
        T_i = κL * κR / abs2((ω[i] - ω_ph + im * 0.5 * κ + χ[i]))
        T[i] = T_i
    end
    return T
end
function linear_transmission(ω, χ, ω_ph, κL, κR, κ)
    κL * κR / abs2((ω - ω_ph + im * 0.5 * κ + χ))
end

function linear_absorption(ω, ω_ph, ω_v, κL, κ, γ, N, g)
    κL * γ * N * g^2 / abs2((ω - ω_ph + im * 0.5 * κ) * (ω - ω_v + im * 0.5 * γ) - N * g^2)
end
