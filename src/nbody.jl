using Plots, DataFrames, CSV, LinearAlgebra
include("utility.jl")

const G = 100.0
const S = 0.05

mutable struct Particle
    m::Float64
    r::Vector{Float64}
    v::Vector{Float64}
end

mutable struct System
    particles::Vector{Particle}
    #cutoff_distance::Float64
    #a::Float64
    #b::Float64
    #c::Float64
    #alpha::Float64
    #beta::Float64
    #gamma::Float64
end

# returns acceleration of p1 due to the force between p1 and p2
function a_g(p1::Particle, p2::Particle)
    r = p2.r - p1.r
    return G * p2.m * r / (norm(r)^2 + S^2)^1.5
end

# returns the net acceleration of a particle due to the forces from all other
# particles within the cutoff distance
function net_a_g(p1::Particle, s::System)
    a = [0.0, 0.0, 0.0]

    for p2 in s.particles
        if p1 != p2
            a += a_g(p1, p2)
        end
    end

    return a
end

# move system s forward by dt using velocity-verlet method
function step(s::System, dt::Float64)
    for p in s.particles
        a1 = net_a_g(p, s)
        p.r += p.v * dt + 0.5 * a1 * dt^2
        a2 = net_a_g(p, s)
        p.v += 0.5 * (a1 + a2) * dt
    end
end

# numerically integrates system particle positions from 0 to t
function integrate(s::System, t::Float64, dt::Float64, save_interval::Float64)
    rows = Int(t / save_interval)
    columns = 3 * length(s.particles) + 1
    data = Array{Float64}(undef, rows, columns)
    data[1, :] = vcat([0.0], [s.particles[j].r[k] for
        j = 1:length(s.particles) for k = 1:3])
    save_counter = 0.0
    i = 2

    for t_i in 0:dt:t
        step(s, dt)
        
        if save_counter >= save_interval
            data[i, :] = vcat([t_i], [s.particles[j].r[k] for
                j = 1:length(s.particles) for k = 1:3])
            i += 1
            save_counter = 0.0
        else
            save_counter += dt
        end
    end

    data = data[1:i - 1, :]
    return data
end

function main()
    system, t, dt, save_interval = load("input.txt")
    data = integrate(system, t, dt, save_interval)
    save(data, "output")
    create_plot(data, "plot")
end

main()