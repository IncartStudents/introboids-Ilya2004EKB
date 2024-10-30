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
    accel::Vector{Tuple{Float64, Float64}}
    height::Float64
    width::Float64
    function WorldState(n_boids, h, w)
        max_vel=1
        vel_module= [(rand() * max_vel ) for _ in 1:n_boids]
        boids = [(rand(0:w), rand(0:h)) for _ in 1:n_boids]
        vel_angle = [(rand() * 2π - π) for _ in 1:n_boids]
        vel_vector = [(vel_module[k] * cos(vel_angle[k]), vel_module[k] * sin(vel_angle[k])) for k in 1:n_boids]
        accel = [(0.0, 0.0) for _ in 1:n_boids]
        new(boids, vel_vector, vel_angle, vel_module, accel, h, w)
    end
end

function bros_detector(position, distance, k, m)
    x0=position[k][1]
    y0=position[k][2]
    x=position[m][1]
    y=position[m][2]
    if (abs(x-x0) ≤ distance) && (abs(y-y0) ≤ distance) && m!=k 
        return [x,y]
    end
end

function velocity(state::WorldState, k)
    state.vel_vector[k] = (state.vel_vector[k][1] - state.accel[k][1], state.vel_vector[k][2] - state.accel[k][2])
    return nothing
end

function cohesion(state::WorldState, position, grp_x, grp_y, k, max_accel)
    sum_x = sum(grp_x)
    sum_y = sum(grp_y)
    avg_x = sum_x / length(grp_x)
    avg_y = sum_y / length(grp_y)
    accel_angle = atan((avg_y - position[k][2])/(avg_x - position[k][1]))
    accel = sqrt((avg_y - position[k][2])^2 + (avg_x - position[k][1])^2) * 0.4
    if accel > max_accel
        accel = max_accel
    end
    return (accel * cos(accel_angle), accel * sin(accel_angle))
end


function borders(state::WorldState, k)
    if state.boids[k][1] ≥ state.width || state.boids[k][1] ≤ 0
        state.vel_vector[k] = (-state.vel_vector[k][1], state.vel_vector[k][2])
    end
    if state.boids[k][2] ≥ state.height || state.boids[k][2] ≤ 0
        state.vel_vector[k] = (state.vel_vector[k][1], -state.vel_vector[k][2])
    end
    return nothing
end

function update!(state::WorldState, n_boids)
    for k in 1:n_boids
        state.boids[k] = state.boids[k] .+ state.vel_vector[k] 
        borders(state, k)
    end

    distance = 10
    weight_coh = 1
    max_accel = 0.05
    for k in 1:n_boids
        groups_x = Float64[]
        groups_y = Float64[]
        for m in 1:n_boids
            bros = bros_detector(state.boids, distance, k, m) 
            if bros !== nothing
                push!(groups_x, bros[1])
                push!(groups_y, bros[2])
            end
        end
        if isempty(groups_x) == 0
            state.accel[k] = (cohesion(state, state.boids, groups_x, groups_y, k, max_accel) .* weight_coh)
            velocity(state, k)
        end
    end
    return nothing
end

function (@main)(ARGS)
    h = 30
    w = 30
    n_boids = 2
    state = WorldState(n_boids, h, w)
    anim = @animate for time = 1:100
        update!(state, n_boids)
        boids = state.boids
        scatter(boids, xlim = (0, state.width), ylim = (0, state.height))
    end
    gif(anim, "boids.gif", fps =10)
end

end
#dew
using .Boids
Boids.main("")
