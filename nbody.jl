using LinearAlgebra
using DataFrames
using CSV
using Plots

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
    a = [0.0, 0.0]

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
function integrate(s::System, t_end::Float64, dt::Float64 = 0.0001,
    save_interval::Float64 = 0.001)

    rows = Int(t_end / save_interval)
    columns = 2 * length(s.particles) + 1
    data = Array{Float64}(undef, rows, columns)
    data[1, :] = vcat([0.0], [p.r[i] for i = 1:2 for p in s.particles])
    save_counter = 0.0
    i = 2

    for t in 0:dt:t_end
        step(s, dt)
        
        if save_counter >= save_interval
            data[i, :] = vcat([t], [p.r[j] for j = 1:2 for p in s.particles])
            i += 1
            save_counter = 0.0
        else
            save_counter += dt
        end
    end

    data = data[1:i - 1, :]
    return data
end

function save(data::Matrix, file_name::String)
    file_to_save = open(file_name * ".csv", "w")
    CSV.write(file_to_save, DataFrame(data, :auto), delim = ",")
    close(file_to_save)
end

function load(file_name::String)
    file = open(file_name, "r")
    system = 1
    close(file)
    return system
end

function create_plot(data::Matrix, file_name::String)
    x = 2
    y = 3
    new_plot = plot(data[:, x], data[:, y])

    for i = 1:((size(data, 2) - 1) / 2) - 1
        x += 2
        y += 2
        plot!(new_plot, data[:, x], data[:, y])
    end

    savefig(new_plot, file_name * ".svg")
end

function main()
    #p1 = Particle(1.0, [0.0, -5.0], [1.5, 0.0])
    #p2 = Particle(10.0, [-2.0, 5.0], [0.0, 0.0])
    #p3 = Particle(3.0, [0.0, 0.0], [-2.0, 0.0])
    #system = System([p1, p2, p3])
    load("input.txt")
    data = integrate(system, 1.2)
    create_plot(data, "output")
    save(data, "output")
end

main()