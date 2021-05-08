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