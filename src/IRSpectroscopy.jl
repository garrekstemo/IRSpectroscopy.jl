module IRSpectroscopy

include("input-output.jl")
include("molecular_potential.jl")

export transmission,
       pp_spectrum,
       morse_potential,
       harmonic,
       energy_levels

end