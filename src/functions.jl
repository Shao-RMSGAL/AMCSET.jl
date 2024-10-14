include("constants.jl")

using ForwardDiff
using Plots
using Roots
using LinearAlgebra
using Distributed

# Reduced mass
#   M₁ = Incident mass (u)
#   M₂ = Target mass (u)
#   return = Reduced mass (u)
Mc(M₁::Real, M₂::Real)::Real = M₁ * M₂ / (M₁ + M₂)

# Center-of-mass (CM) total kinetic energy
#   Mc = Reduced mass (u)
#   V₀ = Incident velocity (m/s)
#   return = Center-of-mass energy (eV)
function Ec(Mc::Real, V₀::Real)::Real
    factor = ustrip(uconvert(u"eV", 1u"u * (m/s)^2"))
    return factor * (1 / 2) * Mc * V₀^2
end

# Incident velocity
# E = Energy (eV)
# M = mass (u)
# return = Velocity (m/s)
function V₀(E::Real, M::Real)::Real
    factor = ustrip(uconvert(u"m/s", 1u"√(eV / u)"))
    return factor * √(2 * E / M)
end

# Free-flying path length, Eq. 7-25 and L = ∛(1/N)
# M₁ = Mass of incident particle (u)
# M₂ = Mass of target particle (u)
# ε = Reduced energy (unitless)
# a = Screening length (Å)
# N = Atom density (atoms/Å³)
# impuslse_valid = Boolean for whether to do impulse approx
# return = Free flight path length (m)
function L(M₁::Real, M₂::Real, ε::Real, a::Real, N::Real, impulse_valid::Bool)::Real
    if ε > 10 && impulse_valid
        E = ε * a
        T_min = displacement_threshold
        b_max_val = b_max(E, T_min, a, M₁, M₂)
        L = 1/((b_max_val * a)^2 * π * N)
    else
        L = ∛(1 / N)
    end
    return L
end

# Impact parameter, Eq. 7-23 and 7-27
# N = atom density (atoms/Å³)
# L = length (Å)
# ε = reduced energy (dimensionless)
# return = Impact parameter (Å)
function p(N::Real, L::Real, ε::Real)::Real
    if ε > 10
        p = abs(-log(rand()) / (π * N * L)) # Eq. 7-23
    else
        p = √(rand() / (π * N^(2 / 3))) # Eq. 7-27
    end
    return p
end

# Maximum reduced impact parameter Eq. 7-37
# E = Reduced energy (eV)
# Eₘᵢₙ = Minimum energy (eV)
# a = Screening length (Å)
# M₁ = Incident mass (u)
# M₂ = Target mass (u)
# return = reduced impact parameter (dimensionless)
function b_max(E::Real, Eₘᵢₙ::Real, a::Real, M₁::Real, M₂::Real)::Real
    γ = 4M₁ * M₂ / (M₁ + M₂)^2
    ε = E / a
    εₘᵢₙ = Eₘᵢₙ / a
    ξ = √(ε * εₘᵢₙ / γ)
    return 1 / (ξ + √ξ + 0.125ξ^0.1)
end

# Screening length (TF)
#   Z₁ = Incident atomic number (no units)
#   Z₂ = Target atomic number (no units)
a(Z₁, Z₂) = (1 / 4) * ∛(9π^2 / 2) * a₀ / (Z₁^(2 / 3) + Z₂^(2 / 3))

# Universal screening function, Eq. 2-74
#   x = Reduced radius (r/aᵤ)
#   return = Screening potential (dimensionsless)
Φᵤ(x::Real)::Real =
    0.1818exp(-3.2x) + 0.5099exp(-0.9423x) + 0.2802exp(-0.4028x) + 0.2817exp(-0.2016x)

# Screening length (Universal). Eq. 2-73
#   Z₁ = Incident atomic number (no units)
#   Z₂ = Target atomic number (no units)
#   return = Screening length (Å)
aᵤ(Z₁::Real, Z₂::Real)::Real = (1 / 4) * ∛(9π^2 / 2) * a₀ / (Z₁^(0.23) + Z₂^(0.23))

# Potential from screening function, Eq. 2-69
# Φ = Screening function as a function of reduced radius(dimensionless) (dimensionless)
# Z = Atomic number (dimensionless)
# r = Separation (Å)
# return = Potential as a function of r(Å) (eV)
function V(Φ₁::Function, Z₁::Real, Z₂::Real)::Function
    factor = ustrip(uconvert(u"eV", 1u"(cm^(3/2) * g^(1/2) / s)^2/Å"))
    return (r) -> factor * Φ₁(r / aᵤ(Z₁, Z₂)) * Z₁ * Z₂ * e_stat^2 / r
end

# Derivative of function V
# V = Potential function (eV)
# r = Radial distance (Å)
# return = Derivative of potential as a function of r(Å) (eV)
function dV(V::Function)::Function
    return (r) -> ForwardDiff.derivative(V, r)
end

# Collision diameter, Eq. 2-31
# Z₁ = Incident Z number (unitless)
# Z₂ = Target Z number (unitless)
# e = Elementary charge (statC)
# Mc = Reduced mass (u)
# V₀ = Incident velocity (m/s)
# return = Collision diameter (Å)
function d(Z₁::Real, Z₂::Real, Mc::Real, V₀::Real)::Real
    factor = ustrip(uconvert(u"Å", 1u"(cm^(3/2) * g^(1/2) / s)^2 / (u * (m/s)^2)"))
    return factor * 4 * Z₁ * Z₂ * e_stat^2 / (Mc * V₀^2)
end

# Closest approach, Eq. 7-2
#   V = Potential function of r(Å) (eV)
#   Ec = Center of mass energy (eV)
#   p = Impact parameter (Å)
#   d = Collision diameter (Å)
#   return = Closest approach (Å)
function r₀(V::Function, Ec::Float64, p::Float64, d::Float64)::Real
    F(r₀) = 1 - V(r₀) / Ec - (p / r₀)^2
    dF(r₀) = ForwardDiff.derivative(F, float(r₀))
    r₀_guess = d / 10
    result = find_zero((F, dF), r₀_guess, Roots.Newton())
    return result
end

# Radius of curvature, Eq. 7-4
# V = Potential as a function of r(Å) (eV)
# dV = Deriavtive of V with respect to r(Å) (eV/Å)
# Ec = Center of mass energy (eV)
# r₀ = Closest approach (Å)
# return = Radius of curvature (Å)
function ρ(V::Function, dV::Function, Ec::Real, r₀::Real)::Real
    return 2abs(Ec - V(r₀)) / (-dV(r₀))
end

# Reduced energy, Eq. 2-83
#   a = Screening length (Å)
#   Ec = Center of mass energy (eV)
#   Z₁ = Incident atomic number (unitless)
#   Z₂ = Target atomic number (unitless)
#   return = Reduced energy (dimensionless)
function ε(a::Real, Ec::Real, Z₁::Real, Z₂::Real)::Real
    result = a * Ec / (Z₁ * Z₂ * e_stat^2)
    factor = ustrip(uconvert(NoUnits, 1u"Å * eV / (cm^(3/2) * g^(1/2) / s)^2"))
    return factor * result
end

# Final scattering angle (Magic formula)
#   p = impact parameter (Å)
#   r₀ = closest approach (Å)
#   ρ = CM ρ₁ + ρ₂ radii of curvature (Å)
#   Δ = correction factor (dimensionless)
#   a = screening distance (Å)
#   return = Scattering angle (radians)
function Θ(p::Real, r₀::Real, ρ::Real, Δ::Real, a::Real)::Real
    B = p / a
    R₀ = r₀ / a
    Rc = ρ / a
    return acos((B + Rc + Δ) / (R₀ + Rc))
end

# Magic formula correction parameter:
#   ε = Reduced energy (eV)
#   p = impact parameter (Å)
#   a = screening length (Å)
#   return = Δ (dimensionless)
function Δ(ε::Real, p::Real, a::Real, r₀::Real)::Real
    C₁, C₂, C₃, C₄, C₅ = C_Δ
    α = 1 + C₁ * ε^(-1 / 2)
    β = (C₂ + √ε) / (C₃ + √ε)
    γ = (C₄ + ε) / (C₅ + ε)
    B = p / a
    A = 2α * ε * B^β
    G = γ / (√(1 + A^2) - A)
    R₀ = r₀ / a
    return A * (R₀ - B) / (1 + G)
end

# High energy sin²(Θ/2) (Only use if ε > 10)
#   ε = reduced enery
#   b = reduced impact parameter
#   return = sin²(Θ/2) where Θ is scattering angle (dimensionsless)
function sin2Θ_2(ε::Real, b::Real)::Real
    return 1 / (1 + (1 + b * (1 + b)) * (2 * ε * b)^2)
end

# Efficient nuclear energy loss
#   M₁ = Incident atom mass (u)
#   M₂ = Target atom mass (u)
#   ε = reduced enery
#   b = reduced impact parameter
#   return = Transferred kinetic energy (eV)
function T_eff(M₁::Real, M₂::Real, E::Real, ε::Real, b::Real)::Real
    return 4 * M₁ * M₂ * E * sin2Θ_2(ε, b) / (M₁ + M₂)^2
end

# Nuclear energy loss
#   M₁ = Incident atom mass (u)
#   M₂ = Target atom mass (u)
#   E = Incident energy (eV)
#   Θ = Scattering angle (rad)
#   return = Transferred kinetic energy (eV)
function T(M₁::Real, M₂::Real, E::Real, Θ::Real)::Real
    return 4 * M₁ * M₂ / (M₁ + M₂)^2 * E * sin(Θ)^2
end

# Laboratory scattering angle
#   Θ = CM Scattering angle (rad)
#   M₁ = Incident mass (u)
#   M₂ = Target mass (u)
#   return = Relative laboratory scattering angle (rad)
ϑ(Θ::Real, M₁::Real, M₂::Real) = atan(sin(Θ) / (cos(Θ) + (M₁ / M₂)))

# Azimuthal scattering angle
# return = Azimuthal scattering angle (rad)
ϕ()::Real = 2π * rand()

# Angle with respect to target normal axis
#   αᵢ₋₁ = Previous normal axis angle (rad)
#   ϑᵢ = Laboratory scattering angle (rad)
#   ϕ = Azimuthal scattering angle (rad)
#   return = Absolute lab angle (rad)
αᵢ(αᵢ₋₁::Real, ϑᵢ::Real, ϕᵢ::Real)::Real =
    acos(cos(αᵢ₋₁) * cos(ϕᵢ) + sin(αᵢ₋₁) * sin(ϑᵢ) * cos(ϕᵢ))

# Convert relative angle to absolute angle for incident and target particle
# α₀ = Initial angle from z-axis (rad)
# ϕ₀ = Initial azimuthal angle (rad)
# ϑ₁ = Deflection angle from velocity vector for incident particle (rad)
# ϑ₁ = Deflection angle from velocity vector for target particle (rad)
# return (θ₁, α₁, θ₂, α₂) New incident particle (1) and target(2) angles (rad)
function relative_to_absolute(α₀, ϕ₀, ϑ₁, ϑ₂)
    α₁ᵣ = rand() * 2π
    α₂ᵣ = α₁ᵣ + π
    θ₀ = α₀
    α₀ = ϕ₀
    θ₁ᵣ = ϑ₁
    θ₂ᵣ = ϑ₂


    X =
        sin(θ₁ᵣ) * cos(α₁ᵣ) * sin(α₀) +
        sin(θ₁ᵣ) * sin(α₁ᵣ) * cos(θ₀) +
        cos(θ₁ᵣ) * sin(θ₀) * cos(α₀)
    Y =
        -sin(θ₁ᵣ) * cos(α₁ᵣ) * cos(α₀) +
        sin(θ₁ᵣ) * sin(α₁ᵣ) * cos(θ₀) +
        cos(θ₁ᵣ) * sin(θ₀) * sin(α₀)
    Z = -sin(θ₁ᵣ) * sin(α₁ᵣ) * sin(θ₀) + cos(θ₁ᵣ) * cos(θ₀)
    if Z > 0
        θ₁ = atan(√(X^2 + Y^2) / Z)
    elseif Z == 0
        θ₁ = π / 2
    else
        θ₁ = π + atan(√(X^2 + Y^2) / Z)
    end

    if sin(θ₁) ≠ 0
        if X > 0
            α₁ = atan(Y / X)
        elseif X == 0
            α₁ = π - (Y > 0 ? 1 : -1) * π / 2
        else
            α₁ = π + atan(Y / X)
        end
    else
        α₁ = 0
    end

    Z = -sin(θ₂ᵣ) * sin(α₂ᵣ) * sin(θ₀) + cos(θ₂ᵣ) * cos(θ₀)
    X =
        sin(θ₂ᵣ) * cos(α₂ᵣ) * sin(α₀) +
        sin(θ₂ᵣ) * sin(α₂ᵣ) * cos(θ₀) +
        cos(θ₂ᵣ) * sin(θ₀) * cos(α₀)
    Y =
        -sin(θ₂ᵣ) * cos(α₂ᵣ) * cos(α₀) +
        sin(θ₂ᵣ) * sin(α₂ᵣ) * cos(θ₀) +
        cos(θ₂ᵣ) * sin(θ₀) * sin(α₀)

    if Z > 0
        θ₂ = atan(√(X^2 + Y^2) / Z)
    elseif Z == 0
        θ₂ = π / 2
    else
        θ₂ = π + atan(√(X^2 + Y^2) / Z)
    end

    if sin(θ₂) ≠ 0
        if X > 0
            α₂ = atan(Y / X)
        elseif X == 0
            α₂ = π - (Y > 0 ? 1 : -1) * π / 2
        else
            α₂ = π + atan(Y / X)
        end
    else
        α₂ = 0
    end
    return θ₁, α₁, θ₂, α₂
end

# Run a simulation
# Z₁ = Incident charge (unitless)
# Z₂ = Target charge (unitless)
# M₁ = Incident mass (u)
# M₂ = Target mass (u)
# E = Energy (eV)
# num = Number of simulations
# return = Vector of vector of positions, one for each particle
    function run_simulation(Z₁, Z₂, M₁, M₂, E_init, ρ_sub, num)
        res_lock = ReentrantLock()
        N = ustrip(AvogadroConstant) * ρ_sub / M₂ / 1E24 # atoms/Å³

        res = Vector{Vector{Tuple{Float64, Float64, Float64}}}()
        s_vec = Vector{Future}()
        V_func = V(Φᵤ, Z₁, Z₂)
        dV_func = dV(V_func)
        a_val = aᵤ(Z₁, Z₂)
        T_min = displacement_threshold
        impulse_valid = true
            
        Mc_val = Mc(M₁, M₂)
        for _ in 1:num
            s = @spawnat :any begin
                pos = (0.0, 0.0, 0.0)
                αᵢ₋₁ = 0
                ϕᵢ₋₁ = 0
                vec = [pos]
                E = E_init
                while (E > displacement_threshold)
                    V₀_val = V₀(E, M₁)
                    Ec_val = Ec(Mc_val, V₀_val)
                    ε_val = ε(a_val, Ec_val, Z₁, Z₂)
                    L_val = L(M₁, M₂, ε_val, a_val, N, impulse_valid)
                    p_val = p(N, L_val, ε_val)
                    b_val = p_val / a_val
                    if ε_val > 10
                        T_val = T_eff(M₁, M₂, E, ε_val, b_val)
                        Θ_val = √(asin(sin2Θ_2(ε_val, b_val))*2)
                    elseif ε_val > T_min
                        d_val = d(Z₁, Z₂, Mc_val, V₀_val)
                        r₀_val = r₀(V_func, Ec_val, p_val, d_val)
                        ρ_val = ρ(V_func, dV_func, Ec_val, r₀_val)
                        Δ_val = Δ(ε_val, p_val, a_val, r₀_val)
                        Θ_val = Θ(p_val, r₀_val, ρ_val, Δ_val, a_val)
                        T_val = T(M₁, M₂, E, Θ_val)
                    else
                    end
                    ϑ_val = ϑ(Θ_val, M₁, M₂)
                    (αᵢ₋₁, ϕᵢ₋₁, αₜ, ϕₜ) = relative_to_absolute(αᵢ₋₁, ϕᵢ₋₁, ϑ_val, 0)
                    if T_val > displacement_threshold
                        println("Knockon created! $T_val: $αₜ, $ϕₜ")
                    end
                    E = E - T_val
                    pos = (
                    pos[1] + L_val * sin(αᵢ₋₁) * cos(ϕᵢ₋₁),
                    pos[2] + L_val * sin(αᵢ₋₁) * sin(ϕᵢ₋₁),
                    pos[3] + L_val * cos(αᵢ₋₁),
                    )
                        push!(vec, pos)
                end
                @lock res_lock push!(res, vec)
            end
            push!(s_vec, s)
        end

        for future in s_vec
            wait(future) 
        end

        return res
    end
