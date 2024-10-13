using PhysicalConstants.CODATA2018
using Unitful

const a₀ = uconvert(
    u"Å",
    ReducedPlanckConstant / (ElectronMass * SpeedOfLightInVacuum * FineStructureConstant),
)

ε_unit = 1u"Å * eV / C^2"

const C_Δ = [
    0.99229  # * √ε_unit,
    0.011615 # * √ε_unit,
    0.007122 # * √ε_unit,
    9.3066   # * ε_unit,
    14.813   # * ε_unit
    ]

const e_stat = ElementaryCharge * 2997924580u"cm^(3/2) * g^(1/2) * s^(-1)/C"

displacement_threshold = 40
