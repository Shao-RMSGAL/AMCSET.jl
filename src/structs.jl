using Gtk4
using Unitful
using PhysicalConstants.CODATA2018
#  using Base.Threads

struct EntryGroup
    text_label::GtkLabel
    entry::GtkEntry
    unit_label::GtkLabel
    default_value::Float64
end

struct Velocity
    energy::Float64
    θ::Float64
    ϕ::Float64
end

mutable struct Substrate
    Z::Float64
    M::Float64
    N::Float64
    Substrate(Z, M, ρ) = new(Z, M, ρ * ustrip(AvogadroConstant) / M / 1e24) # atoms/Å³
end

const Coordinate = Tuple{Float64, Float64, Float64}

mutable struct Particle
    position::Coordinate
    velocity::Velocity
    Z::Float64
    M::Float64
    Particle(position, velocity, Z, M) = new(position, velocity, Z, M)
    Particle(position, velocity, substrate::Substrate) =
        new(position, velocity, substrate.Z, substrate.M)
end

#  # Define our thread-safe structure
#  struct ThreadSafeVectorOfVectors
#      data::Vector{Vector{Coordinate}}
#      lock::ReentrantLock
#      counter::Atomic{Int}
#  end

#  # Constructor
#  function ThreadSafeVectorOfVectors()
#      ThreadSafeVectorOfVectors(Vector{Vector{Coordinate}}(), ReentrantLock(), Atomic{Int}(0))
#  end
