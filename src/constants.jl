using PhysicalConstants.CODATA2018
using Unitful

# Bohr radius in Å
const a₀ = ustrip(uconvert(
    u"Å",
    ReducedPlanckConstant / (ElectronMass * SpeedOfLightInVacuum * FineStructureConstant)
    ))


const C_Δ = [
    0.99229  
    0.011615 
    0.007122 
    9.3066   
    14.813   
    ]

# Charge of an electron in statC
const e_stat = ustrip(ElementaryCharge * 2997924580u"cm^(3/2) * g^(1/2) * s^(-1)/C")

displacement_threshold = 40
