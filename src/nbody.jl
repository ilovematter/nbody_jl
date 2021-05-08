using Plots, DataFrames, CSV, LinearAlgebra

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

function save(data::Matrix, file_name::String)
    file_to_save = open(file_name * ".csv", "w")
    CSV.write(file_to_save, DataFrame(data, :auto), delim = ",")
    close(file_to_save)
end

function load(file_name::String)
    file = open(file_name, "r")
    system = System([])
    t = 0.0
    dt = 0.0
    save_interval = 0.0

    for line in readlines(file)
        line = split(line, " ")

        if cmp(line[1], "time") == 0
            t = parse(Float64, line[2])
        elseif cmp(line[1], "step_size") == 0
            dt = parse(Float64, line[2])
        elseif cmp(line[1], "save_interval") == 0
            save_interval = parse(Float64, line[2])
        elseif cmp(line[1], "particle") == 0
            m = parse(Float64, line[2])
            r = [parse(Float64, line[3]), parse(Float64, line[4]),
                parse(Float64, line[5])]
            v = [parse(Float64, line[6]), parse(Float64, line[7]),
                parse(Float64, line[8])]
            push!(system.particles, Particle(m, r, v))
        else
            throw("syntax error")
        end
    end

    close(file)
    return system, t, dt, save_interval
end

function create_plot(data::Matrix, file_name::String)
    x = 2
    y = 3
    z = 4
    new_plot = plot(data[:, x], data[:, y], data[:, z], legend = false,
        aspect_ratio = :equal)

    for i = 1:((size(data, 2) - 1) / 3) - 1
        x += 3
        y += 3
        z += 3
        plot!(new_plot, data[:, x], data[:, y], data[:, z], legend = false,
            aspect_ratio = :equal)
    end

    savefig(new_plot, file_name * ".svg")
end

function main()
    system, t, dt, save_interval = load("input.txt")
    data = integrate(system, t, dt, save_interval)
    save(data, "output")
    create_plot(data, "plot")
end

main()