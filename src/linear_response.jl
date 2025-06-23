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
# function linear_transmission(ω, ω_ph, ω_v, N, g, κL, κR, κ, γ)
#     κL * κR * ((ω - ω_v)^2 + (0.5 * γ)^2) / abs2((ω - ω_ph + im * 0.5 * κ) * (ω - ω_v + im * 0.5 * γ) - N * g^2)
# end
function linear_transmission(ω, ω_ph, ω_v, N, g, κL, κR, κ, γ)
    κL * κR / abs2((ω - ω_ph + im * 0.5 * κ + χ(ω, ω_v, N, g, γ)))
end

function linear_transmission_anharmonic(ω, ω_ph, ω_v1, ω_v2, N, g, transfer, κL, κR, κ, γ1, γ2)
    N1 = N * (1 - transfer)
    N2 = N * transfer
    κL * κR / abs2((ω - ω_ph + im * 0.5 * κ + χ(ω, ω_v1, N1, g, γ1) + χ(ω, ω_v2, N2, g, γ2)))
end

function linear_absorption(ω, ω_ph, ω_v, κL, κ, γ, N, g)
    κL * γ * N * g^2 / abs2((ω - ω_ph + im * 0.5 * κ) * (ω - ω_v + im * 0.5 * γ) - N * g^2)
end
