# function χ(ω, ω_v, N, g, γ)
#     -N * abs2(g) / (ω - ω_v + im * 0.5 * γ)
# end
function χ(ω, p::Vector, μ, ω_zy, N, γ, V_mol)
    a = zeros(ComplexF64, length(ω))
    for i in eachindex(ω)
        for y in 1:length(p)
            for z in 1:length(p)
                a[i] += (p[y] - p[z]) * abs2(μ[z, y]) / (ω[i] - ω_zy[z, y] + im * 0.5 * γ)
            end
        end
    end
    return a .* -N / V_mol
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
function linear_transmission(ω, χ::Vector, ω_ph, κL, κR, κ)
    Ts = similar(ω, Float64)
    for i in eachindex(ω)
        Ts[i] = κL * κR / abs2((ω[i] - ω_ph + im * 0.5 * κ + χ[i]))
    end
    return Ts
end
function linear_transmission(ω, χ, ω_ph, κL, κR, κ)
    κL * κR / abs2((ω - ω_ph + im * 0.5 * κ + χ))
end

function linear_absorption(ω, ω_ph, ω_v, κL, κ, γ, N, g)
    κL * γ * N * g^2 / abs2((ω - ω_ph + im * 0.5 * κ) * (ω - ω_v + im * 0.5 * γ) - N * g^2)
end
