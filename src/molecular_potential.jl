
function morse_potential(q, q_eq, α, De)
    @. De * (1 - exp(-α * (q - q_eq)))^2
end

function harmonic(q, q_eq, α, De)
    @. De * (-α * (q - q_eq))^2
end

function energy_levels(n, ν_0, ν_exe)
    return @. (n + 0.5) * ν_0 - ν_exe * (n + 0.5)^2
end