include("constants.jl")

# Reduced mass
#   M₁ = Incident mass (mass)
#   M₂ = Target mass (mass)
Mc(M₁, M₂) = M₁ * M₂ / (M₁   + M₂)

# Center-of-mass (CM) total kinetic energy 
#   Mc = Reduced mass (mass)
#   V₀ = Incident velocity (speed) 
Ec(Mc, V₀) = (1/2) * Mc * V₀^2

# Center-of-mass totla kinetic energy alternative
#   E = Incident energy (energy)
#   M₁ = Incident mass (mass)
#   M₂ = Target mass (mass)
Ec(E, M₁, M₂) = E / (1 + M₁/M₂)

# Screening length (TF)
#   Z₁ = Incident atomic number (no units)
#   Z₂ = Target atomic number (no units)
a(Z₁, Z₂) = (1/4) * ∛(9π^2/2) * a₀ / (Z₁^(2/3) + Z₂^(2/3))

# Reduced energy
#   a = Screening length (distance)
#   Ec = Center of mass energy (energy)
#   Z₁ = Incident atomic number (no units)
#   Z₂ = Target atomic number (no units)
ε(a, Ec, Z₁, Z₂) = upreferred(a * Ec / (Z₁ * Z₂ * e_stat ^ 2))

# Final scattering angle (Magic formula)
#   p = impact parameter (distance)
#   r₀ = closest approach (distance)
#   ρ = CM ρ₁ + ρ₂ radii of curvature (distance)
#   δ = correction factor (distance)
function Θ(p, r₀, ρ, Δ, a)
    B = p / a
    R₀ = r₀ / a
    Rc = ρ / a
    return acos((B + Rc + Δ)/(R₀ + Rc)) * 2u"rad"
end

# Magic formula correction parameter:
#   ε = Reduced energy (energy)
#   p = impact parameter (distance)
#   a = screening length (distance)
function Δ(ε, p, a, r₀)
    C₁, C₂, C₃, C₄, C₅ = C_Δ
    α = 1 + C₁*ε^(-1/2)
    β = (C₂ + √ε) / (C₃ + √ε)
    γ = (C₄ + ε) / (C₅ + ε)
    B = p / a
    A = 2α*ε*B^β
    G = γ/(√(1 * A^2) - A)
    R₀ = r₀ / a
    return A * (R₀ - B) / (1 + G)
end

# High energy sin²(Θ/2) (Only use if ε > 10)
#   ε = reduced enery
#   b = reduced impact parameter
function sin2Θ_2(ε, b)
    return 1 / (1 + (1 + b * (1 + b)) * (2 * ε * b)^2)
end

# Nuclear energy loss
#   M₁ = Incident atom mass (mass)
#   M₂ = Target atom mass (mass)
#   Θ = Scattering angle (angle)
T(M₁, M₂, E, Θ) = 4 * M₁ * M₂ / (M₁ + M₂)^2 * E * sin(Θ)^2

# Laboratory scattering angle
#   Θ = CM Scattering angle (angle)
#   M₁ = Incident mass (mass)
#   M₂ = Target mass (mass)
ϑ(Θ, M₁, M₂) = atan(sin(Θ) / (cos(Θ) + (M₁/M₂)))u"rad"

# Azimuthal scattering angle
ϕ() = 2π*rand()u"rad"

# Angle with respect to target normal axis
#   αᵢ₋₁ = Previous normal axis angle
#   ϑᵢ = Laboratory scattering angle
#   ϕ = Azimuthal scattering angle
αᵢ(αᵢ₋₁, ϑᵢ, ϕᵢ) = acos(cos(αᵢ₋₁) * cos(ϕᵢ) + sin(αᵢ₋₁) * sin(ϑᵢ)*cos(ϕᵢ))u"rad"

# Impact parameter (high energy)
#   N = Number density of atoms (inverse volume)
#   L = Free flight path length
p(N, L) = abs(-log(rand())/ (π * N * L))

# Impact parameter (low energy)
#   N = Number density of atoms (inverse volume)
p(N) = √abs(rand() / (π * N^(2/3)))

# Free flight length
#   M₁ = Incident particle mass number
#   M₂ = Target particle mass number
#   N = Atomic density
#   a = screening length
#
function L(M₁, M₂, N, a, ε)
    val = uconvert(u"Å", (0.02 * (1 + ustrip(uconvert(u"u", M₁ + M₂)))^2 * ε^2 + 0.1 * ε^1.38) /
    (4π * a^2 * N * log(1 + ε)))
end
