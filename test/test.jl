include("../src/functions.jl")

using PhysicalConstants.CODATA2018

function test()
    print_output = false
    Z₁ = 1
    Z₂ = 3
    M₁ = 1 # u
    M₂ = 6 # u
    E = 100 # eV
    ρ_sub = 1 # g/cm³
    αᵢ₋₁ = 0
    ϕᵢ₋₁ = 0
    pos = (0.0, 0.0, 0.0)
    N = ustrip(AvogadroConstant) * ρ_sub / M₂ / 1E24 # atoms/Å³

    print_output ? println("Test parameters:\nZ₁, M₁: $Z₁, $M₁\nZ₂, M₂: $Z₂, $M₂\nN: $N") : 

    vec = [pos]

    while (E > 40)
        V₀_val = V₀(E, M₁)
        print_output ? println("Incident velocity: $V₀_val m/s") : 
        Mc_val = Mc(M₁, M₂)
        print_output ? println("Reduced mass: $Mc_val u") :
        Ec_val = Ec(Mc_val, V₀_val)
        print_output ? println("Reduced energy: $Ec_val eV") :
        V_func = V(Φᵤ, Z₁, Z₂)
        print_output ? println("Potential at 1 Å: $(V_func(1)) eV") :
        dV_func = dV(V_func)
        print_output ? println("Differential potential at 1 Å: $(dV_func(1)) eV/Å") :
        a_val = aᵤ(Z₁, Z₂)
        print_output ? println("Screening length: $a_val Å") :
        ε_val = ε(a_val, Ec_val, Z₁, Z₂)
        print_output ? println("Reduced energy: $ε_val") :
        L_val = L(M₁, M₂, ε_val, a_val, N)
        print_output ? println("Free-flying length: $L_val Å") :
        p_val = p(N, L_val, ε_val)
        print_output ? println("Impact parameter: $p_val Å") :
        d_val = d(Z₁, Z₂, Mc_val, V₀_val)
        print_output ? println("Collision diameter: $d_val Å") :
        r₀_val = r₀(V_func, Ec_val, p_val, d_val)
        print_output ? println("Minimum approach: $r₀_val Å") :
        ρ_val = ρ(V_func, dV_func, Ec_val, r₀_val)
        print_output ? println("Radius of curvature: $ρ_val Å") :
        Δ_val = Δ(ε_val, p_val, a_val, r₀_val)
        print_output ? println("Magic formula correction: $Δ_val") :
        Θ_val = Θ(p_val, r₀_val, ρ_val, Δ_val, a_val)
        print_output ? println("CM Scattering angle Θ: $Θ_val rad") :
        T_val = T(M₁, M₂, E, Θ_val)
        E = E - T_val
        print_output ? println("Energy transferred: $T_val eV\nNew E: $E") :
        ϑ_val = ϑ(Θ_val, M₁, M₂)
        print_output ? println("Lab. Relative Scattering angle: $ϑ_val") :
        (αᵢ₋₁, ϕᵢ₋₁, _, _) = relative_to_absolute(αᵢ₋₁, ϕᵢ₋₁, ϑ_val, 0)
        print_output ? println("Absolute αᵢ: $αᵢ₋₁") :
        print_output ? println("Absolute ϕᵢ: $ϕᵢ₋₁") :
        pos = (
            pos[1] + L_val * sin(αᵢ₋₁) * cos(ϕᵢ₋₁),
            pos[2] + L_val * sin(αᵢ₋₁) * sin(ϕᵢ₋₁),
            pos[3] + L_val * cos(αᵢ₋₁),
        )
        print_output ? println("New position: $pos") :
        push!(vec, pos)
    end

    return vec
end

res = Vector{Vector{Tuple{Float64, Float64, Float64}}}()
for i in (1:1000)
    push!(res, test())
end

plot(legend = false)
for data in res
    plot!(data)
end
gui()
