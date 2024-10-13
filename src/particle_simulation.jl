include("structs.jl")
include("constants.jl")

using SharedArrays
using Distributed

# Calculate electronic energy loss
function electronic_energy_loss(incident_particle::Particle, substrate::Substrate)
    Z₁ = incident_particle.Z
    Z₂ = substrate.Z
    M₁ = incident_particle.M
    M₂ = substrate.M
    E = incident_particle.velocity.energy
    N = substrate.N

    L = √substrate.N
    Kₗ = 1.212 * Z₁^(7 / 6) * Z₂ / ((Z₁^(2 / 3) + Z₂^(2 / 3))^(3 / 2) * √M₁)
    Sₑ = Kₗ * √E
    Eₚ = E - 1.59 * L * N * Sₑ / (1000)
    return Eₚ
end

impact_parameter(substrate::Substrate) = √(rand() / (π * substrate.N^(2 / 3)))

# Perform a bombardment
function bombardment(incident_particle::Particle, substrate::Substrate, shared_array)
    println("Bombardment: $incident_particle")
    energy = incident_particle.velocity.energy
    L = √substrate.N
    coordinates = Vector{Tuple{Float64,Float64,Float64}}()
    while energy > displacement_threshold
        energy = incident_particle.velocity.energy
        θ = incident_particle.velocity.θ
        ϕ = incident_particle.velocity.ϕ
        x, y, z = incident_particle.position
        push!(coordinates, (x, y, z))

        Eₚ = electronic_energy_loss(incident_particle, substrate)
        P = impact_parameter(substrate)
        θ₁ᵣ, θ₂ᵣ, re = recoil(incident_particle, substrate, P)
        θ₁, ϕ₁, θ₂, ϕ₂ = amagic(incident_particle, θ₁ᵣ, θ₂ᵣ)
        incident_particle.velocity = Velocity(Eₚ - re, θ₁, ϕ₁)
        incident_particle.position =
            (x + L * sin(θ) * cos(ϕ), y + L * sin(θ) * sin(ϕ), z + L * cos(θ))
        target_velocity = Velocity(re, θ₂, ϕ₂)
        if re > displacement_threshold
            target_position = incident_particle.position
            target_particle = Particle(target_position, target_velocity, substrate)
            bombardment(target_particle, substrate, shared_array)
        end
    end

    push!(shared_array, coordinates)
    return shared_array
end

function F(X, columbiavk, AU)
    if X == 0
        return 0
    else
        return columbiavk *
               X *
               (
                   0.35 * exp(-0.3 / X / AU) +
                   0.55 * exp(-1.2 / X / AU) +
                   0.1 * exp(-6 / X / AU)
               )
    end
end

function DF(X, columbiavk, AU)
    DF =
        F(X, columbiavk, AU) / X +
        columbiavk / X * (
            0.35 * exp(-0.3 / X / AU) * 0.3 / AU +
            0.55 * exp(-1.2 / X / AU) * 1.2 / AU +
            0.1 * exp(-6 / X / AU) * 6 / AU
        )
end

function recoil(incident_particle::Particle, substrate::Substrate, impact_parameter)

    θ₁ᵣ = 0
    θ₂ᵣ = 0

    Z₁ = incident_particle.Z
    Z₂ = substrate.Z
    M₁ = incident_particle.M
    M₂ = substrate.M
    P = impact_parameter
    E = incident_particle.velocity.energy

    columbiavk = 0.0143992 * Z₁ * Z₂ # Not clear
    MC = M₁ * M₂ / (M₁ + M₂) # Reduced mass
    V = √(E * 2 / M₁) # Laboratory velocity
    EC = (1 / 2) * MC * V^2 # Center-of-mass energy
    AU = 0.8854 * 0.529 / (√Z₁ + √Z₂)^(2 / 3) # Universal screening distance
    elinhard = EC * AU / columbiavk # Not clear

    # Iteration to find minimum
    AA = P^2
    BB = columbiavk / EC
    CC = -1
    columrmin = (1 / 2) / AA * (-BB + √(BB^2 - 4 * AA * CC))
    caltime = 1

    rmintry1 = columrmin # Assume columb minimum
    rmintry2 = 0
    counts = 0
    while abs(rmintry2 - rmintry1) > 1e-2
        rmintry1 = columrmin
        dv = abs(-DF(rmintry1, columbiavk, AU) / EC - 2 * P^2 * rmintry1) # Derivative of some function
        if (abs(dv) < 1e-6)
            dv = 0.1 # Purpose not clear
        end
        rmintry2 = rmintry1 + (1 - F(rmintry1, columbiavk, AU) / EC - P^2 * rmintry1^2) / dv
        columrmin = rmintry1
        counts += 1
        if counts > 100
            break
        end
    end
    rmin = (rmintry1 + rmintry2) / 2

    rbiersack = 2 * (EC - F(rmin, columbiavk, AU)) / rmin^2 / DF(rmin, columbiavk, AU)
    bbiersack = P / AU
    robiersack = 1 / (rmin * AU)
    rcbiersack = rbiersack / AU

    c1biersack = 0.6743
    c2biersack = 0.009611
    c3biersack = 0.005175
    c4biersack = 10.0
    c5biersack = 6.314

    althabiersack = 1 + c1biersack / √elinhard
    beltabiersack = (c2biersack + √elinhard) / (c3biersack + √elinhard)
    gamabiersack = (c4biersack + elinhard) / (c5biersack + elinhard)
    abiersack = 2 * althabiersack * elinhard * bbiersack^beltabiersack
    gbiersack = gamabiersack / (√(1 + abiersack^2) - abiersack)
    deltabiersack = abiersack * (robiersack - bbiersack) / (1 + gbiersack)


    if P == 0
        calpha1 = π
    else
        calpha1 =
            2 * atan(
                √(
                    (robiersack + rcbiersack)^2 /
                    (bbiersack + rcbiersack + deltabiersack)^2 - 1,
                ),
            )
    end

    cosplusmass = cos(calpha1) + M₁ / M₂
    if cosplusmass < 0
        θ₁ᵣ = π / 2
    elseif cosplusmass < 0
        θ₁ᵣ = π + atan(sin(calpha1) / cosplusmass)
    else
        θ₁ᵣ = atan(sin(calpha1) / cosplusmass)
    end

    re = 4 * EC * MC / M₁ * sin(calpha1 / 2) * sin(calpha1 / 2)
    θ₂ᵣ = (π - calpha1) / 2

    return θ₁ᵣ, θ₂ᵣ, re
end

function amagic(incident_particle, θ₁ᵣ, θ₂ᵣ)
    α₁ᵣ = rand() * 2π
    α₂ᵣ = α₁ᵣ + π
    θ₀ = incident_particle.velocity.θ
    α₀ = incident_particle.velocity.ϕ

    X1 = sin(θ₁ᵣ) * cos(α₁ᵣ)
    Y1 = sin(θ₁ᵣ) * sin(α₁ᵣ)
    Z1 = cos(θ₁ᵣ)

    Y0 = Y1 * cos(θ₀) + Z1 * sin(θ₀)
    Z0 = -Y1 * sin(θ₀) + Z1 * cos(θ₀)
    X0 = X1

    Z = Z0
    X = X0 * sin(α₀) + Y0 * cos(α₀)
    Y = -X0 * cos(α₀) + Y0 * sin(α₀)
    if Z > 0
        θ₁ = atan(√(X^2 + Y^2) / Z)
    elseif Z == 0
        θ₁ = π / 2
    else
        θ₁ = π + atan(√(X^2 + Y^2) / Z)
    end

    if sin(θ₁) ≠ 0
        if X > 0
            α₁ = atan(Y / Z)
        elseif X == 0
            α₁ = π - (Y > 0 ? 1 : -1) * π / 2
        else
            α₁ = π + atan(Y / X)
        end
    else
        α₁ = 0
    end

    X1 = sin(θ₂ᵣ) * cos(α₂ᵣ)
    Y1 = sin(θ₂ᵣ) * sin(α₂ᵣ)
    Z1 = cos(θ₂ᵣ)

    Y0 = Y1 * cos(θ₀) + Z1 * sin(θ₀)
    Z0 = -Y1 * sin(θ₀) + Z1 * cos(θ₀)
    X0 = X1

    Z = Z0
    X = X0 * sin(α₀) + Y0 * cos(α₀)
    Y = -X0 * cos(α₀) + Y0 * cos(α₀)

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

# Launch an ion
#   Z₁ = Incident ion Z-number
#   Z₂ = Target ion Z-number
#   M₁ = Incident ion mass-number
#   M₂ = Target ion mass-numer
#   N = Target density
function initialize(Z₁, Z₂, M₁, M₂, E, ρ)
    is_electron = false
    initial_position = (0.0, 0.0, 0.0)
    initial_velocity = Velocity(E, 0.0, 0.0)
    incident_particle = Particle(initial_position, initial_velocity, Z₁, M₁)
    target = Substrate(Z₂, M₂, ρ)

    shared_array = Vector{Vector{Tuple{Float64,Float64,Float64}}}()

    result = bombardment(incident_particle, target, shared_array)
    return result
end

initialize(20, 3, 50, 6, 100, 1)
