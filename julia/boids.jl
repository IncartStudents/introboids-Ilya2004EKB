module Boids
 import Pkg; Pkg.add("Plots")
using Plots

mutable struct WorldState
    boids::Vector{Tuple{Float64, Float64}}
    vel_vector::Vector{Tuple{Float64, Float64}}
    vel_angle::Vector{Float64}
    height::Float64
    width::Float64
    function WorldState(n_boids, h, w)
        vel_module=0.1
        boids = [(rand(0:w), rand(0:h)) for _ in 1:n_boids]
        vel_angle = [(rand() * 2π - π) for _ in 1:n_boids]
        vel_vector = [(vel_module * cos(vel_angle[k]), vel_module * sin(vel_angle[k])) for k in 1:n_boids]
        new(boids, vel_vector, vel_angle, h, w)
    end
end

function cohesion(position, n_boids, distance)
    groups = [Array{Float64}[] for k in 1:n_boids]
    # boids_x=sort(position[1])
    # boids_y=sort(position[2])
    println(typeof(groups))
    println()
    for k in 1:n_boids
        x0=position[k][1]
        y0=position[k][2]
        for m in 1:n_boids
            x=position[m][1]
            y=position[m][2]
            if (abs(x-x0) ≤ distance) && (abs(y-y0) ≤ distance) && m!=k 
                neighbor = [x,y]
                push!(groups, neighbor)
            end
        end
        println(group)
    end
    return groups
end




function update!(state::WorldState, n_boids)
    for k in 1:n_boids
        state.boids[k] = state.boids[k] .+ state.vel_vector[k] 
        cohesion(state.boids, n_boids, 5)
    end
    return nothing
end

function (@main)(ARGS)
    h = 30
    w = 30
    n_boids = 10
    state = WorldState(n_boids, h, w)
    anim = @animate for time = 1:100
        update!(state, n_boids)
        boids = state.boids
        scatter(boids, xlim = (0, state.width), ylim = (0, state.height))
    end
    gif(anim, "boids.gif", fps = 10)
end

end
#dew
using .Boids
Boids.main("")
