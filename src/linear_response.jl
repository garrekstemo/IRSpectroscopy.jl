function χ(ω)
end

function ξ(ω, ω_ph, κ)
   abs2(ω - ω_ph + im * 0.5 * κ + 0.5 * ω_ph * χ(ω))
end

function transmission(ω, ω_ph, κL, κR, κ)
    κL * κR / ξ(ω, ω_ph, κ)
end

function absorption(ω, ω_ph, κL, κ)
    κL * ω_ph * Fχ(ω) / ξ(ω, ω_ph, κ)
end
function reflection(ω, ω_ph, κL, κR, κ) 
    (ξ(ω, ω_ph, κ) - κL * κR - κL * ω_ph * Fχ(ω)) / ξ(ω, ω_ph, κ)
end