function b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, A)

    bout = Vector{ComplexF64}(undef, length(ω))

    for (i, ω_i) in enumerate(ω)
        bout_i = -im * 0.5 * κ * (ω_i - ω_0 + im * γ_1) * (ω_i - ω_12 + im * γ_3) / ( (ω_i - ω_c + im * 0.5 * κ) * (ω_i - ω_0 + im * γ_1) * (ω_i - ω_12 + im * γ_3) - N * A[i] )
        bout[i] = bout_i
        # push!(bout, bout_i) 
    end

    return bout
end

function A(ω, ω_0, g_1, Δ, γ_m, f_pu, B)

    a = Vector{ComplexF64}(undef, length(ω))

    for (i, ω_i) in enumerate(ω)
        a_i = g_1^2 * (ω_i - (ω_0 - 2*Δ) + 3 * 0.5 * im * γ_m) + f_pu * B[i]
        a[i] = a_i
    end
    return a
end

function B(ω, ω_0, g_1, g_3, Δ, γ_m)

    b = Vector{ComplexF64}(undef, length(ω))

    for (i, ω_i) in enumerate(ω)
        b_i = g_3 * (2 * (2 * g_1 + g_3) * (ω_i - ω_0) + (4 * g_1 + g_3) * im * γ_m) - 4 * Δ * g_1^2
        b[i] = b_i
    end
    b
end

function calculate_ΔT(ω, ω_0, ω_c, ω_12, g_1, g_3, γ_m, γ_1, γ_3, κ, Δ, N, f_pu)
    B_ω = B(ω, ω_0, g_1, g_3, Δ, γ_m)
    A_ω = A(ω, ω_0, g_1, Δ, γ_m, f_pu, B_ω)
    T_fpu = abs2.(b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, A_ω))

    A0 = A(ω, ω_0, g_1, Δ, γ_m, 0.0, B_ω)
    T_0 = abs2.(b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, A0))

    return T_fpu .- T_0
end

function calculateT(ω, ω_0, ω_c, ω_12, g_1, g_3, γ_m, γ_1, γ_3, κ, Δ, N, f_pu)
    B_ω = B(ω, ω_0, g_1, g_3, Δ, γ_m)
    A_ω = A(ω, ω_0, g_1, Δ, γ_m, f_pu, B_ω)
    T = abs2.(b_out(ω, ω_0, ω_c, ω_12, γ_1, γ_3, κ, N, A_ω))

    return T
end

# function calculate_ΔT()
# end

