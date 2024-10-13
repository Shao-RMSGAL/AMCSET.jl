include("functions.jl")

function test()
    Z₁ = 1
    M₁ = 1u"u"
    Z₂ = 5
    M₂ = 10u"u"
    E = 1u"MeV"
    N = 8e22u"cm^-3"
    r₀ = 10u"Å" 
    p = 5u"Å"
    ρ = 1u"Å"
    αᵢ₋₁ = 0u"rad"

    a_val = a(Z₁, Z₂)
    println("a: $a_val") 
    Ec_val = Ec(E, M₁, M₂)
    println("Ec val: $Ec_val")
    ε_val = ε(a_val, Ec_val, Z₁, Z₂)
    println("ε val: $ε_val")
    del = Δ(ε_val, p, a_val, r₀)
    println("Delta: $del")
    Θ_val = Θ(p, r₀, ρ, del, a_val)
    println("Magic Θ: $Θ_val")
    T_val = T(M₁, M₂, E, Θ_val)
    println("Energy loss: $T_val")
    ϑ_val = ϑ(Θ_val, M₁, M₂)
    println("ϑ: $ϑ_val")
    ϕ_val = ϕ()
    println("ϕ: $ϕ_val")
    αᵢ_val = αᵢ(αᵢ₋₁, ϑ_val, ϕ_val)
    println("New angle: $αᵢ_val")
    L_val = L(M₁, M₂, N, a_val, ε_val)
    println("L value: $L_val")
end

test()
