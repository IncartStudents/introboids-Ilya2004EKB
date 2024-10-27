module Boids
 import Pkg; Pkg.add("Plots")
 import Statistics
using Plots
using Statistics

mutable struct WorldState
    boids::Vector{Tuple{Float64, Float64}}
    vel_vector::Vector{Tuple{Float64, Float64}}
    vel_angle::Vector{Float64}
    vel_module::Vector{Float64}
    height::Float64
    width::Float64
    function WorldState(n_boids, h, w)
        vel_module= [(rand()*0.7) for _ in 1:n_boids]
        boids = [(rand(0:w), rand(0:h)) for _ in 1:n_boids]
        vel_angle = [(rand() * 2π - π) for _ in 1:n_boids]
        vel_vector = [(vel_module[k] * cos(vel_angle[k]), vel_module[k] * sin(vel_angle[k])) for k in 1:n_boids]
        new(boids, vel_vector, vel_angle, vel_module, h, w)
    end
end

function bros_detector(position, n_boids, distance, k, m)
    x0=position[k][1]
    y0=position[k][2]
    x=position[m][1]
    y=position[m][2]
    if (abs(x-x0) ≤ distance) && (abs(y-y0) ≤ distance) && m!=k 
        return [x,y]
    end
end

function velocity(state::WorldState, adj_angle, k, weight_coh)
    state.vel_angle[k]= state.vel_angle[k] - adj_angle * weight_coh
    state.vel_vector[k] = state.vel_module[k] * cos(state.vel_angle[k]), state.vel_module[k] * sin(state.vel_angle[k])
    return nothing
end

function cohesion(position, grp_x, grp_y, k)
    sum_x = sum_y = 0
    for i in eachindex(grp_x)
        sum_x = sum_x + grp_x[i]
        sum_y = sum_y + grp_y[i]
    end
    avg_x = sum_x / length(grp_x)
    avg_y = sum_y / length(grp_y)
    angle = atan((avg_y - position[k][2])/(avg_x - position[k][1]))
    return angle
end



function update!(state::WorldState, n_boids)
    for k in 1:n_boids
        state.boids[k] = state.boids[k] .+ state.vel_vector[k] 
    end

    distance = 6
    weight_coh = 0.4
    for k in 1:n_boids
        groups_x = []
        groups_y = []
        for m in 1:n_boids
            bros = bros_detector(state.boids, n_boids, distance, k, m) 
            if bros !== nothing
                push!(groups_x, bros[1])
                push!(groups_y, bros[2])
            end
        end
        if isempty(groups_x) == 0
            adj_angle = (cohesion(state.boids, groups_x, groups_y, k) - state.vel_angle[k])
            velocity(state, adj_angle, k, weight_coh)
        end
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
