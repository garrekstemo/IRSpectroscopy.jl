function dielectric_real(ω, A, ω_0, Γ)
    return @. A * (ω_0^2 - ω^2) / ((ω^2 - ω_0^2)^2 + (Γ * ω)^2)
end

function dielectric_imag(ω, A, ω_0, Γ)
    return @. A * Γ * ω / ((ω^2 - ω_0^2)^2 + (Γ * ω)^2)
end

function nk(ω, A, ω_0, Γ, n_bg = 1.0)
    ε1 = dielectric_real.(ω, A, ω_0, Γ) .+ n_bg^2
    ε2 = dielectric_imag.(ω, A, ω_0, Γ)
    n = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) + ε1) / 2)
    k = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) - ε1) / 2)
    return n, k
end

function nk2(ω, p)
    A1, ω1, Γ1, A2, ω2, Γ2, n_bg = p
    ε1 = dielectric_real.(ω, A1, ω1, Γ1) .+ dielectric_real.(ω, A2, ω2, Γ2) .+ n_bg^2
    ε2 = dielectric_imag.(ω, A1, ω1, Γ1) .+ dielectric_imag.(ω, A2, ω2, Γ2)
    n = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) + ε1) / 2)
    k = @. sqrt((sqrt(abs2(ε1) + abs2(ε2)) - ε1) / 2)
    return n, k
end
