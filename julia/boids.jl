module Boids
 import Pkg; Pkg.add("Plots")
 import Statistics
using Plots
using Statistics

mutable struct WorldState
    boids::Vector{Tuple{Float64, Float64}}
    vel_vector::Vector{Tuple{Float64, Float64}}
    accel::Vector{Tuple{Float64, Float64}}
    height::Float64
    width::Float64
    max_vel::Float64
    max_accel::Float64
    distance::Float64
    weight_al::Float64
    weight_coh::Float64
    weight_sep::Float64
    n_boids::Int64
    function WorldState(h, w)
        max_vel = 1 
        n_boids = 2
        max_accel = 1
        distance = 8
        weight_coh = 0
        weight_sep = 0
        weight_al = 0.9
        boids = [(rand(0:w), rand(0:h)) for _ in 1:n_boids]
        vel_angle = [(rand() * 2π - π) for _ in 1:n_boids]
        vel_vector = [(max_vel * cos(vel_angle[k]), max_vel * sin(vel_angle[k])) for k in 1:n_boids]
        accel = [(0.0, 0.0) for _ in 1:n_boids]
        new(boids, vel_vector, accel, h, w, max_vel, max_accel, distance, 
        weight_al, weight_coh, weight_sep,n_boids)
    end
end

function bros_detector(state::WorldState, k, m)
    position = state.boids
    x0=position[k][1]
    y0=position[k][2]
    x=position[m][1]
    y=position[m][2]
    if (sqrt((x - x0)^2 + (y - y0)^2) ≤ state.distance) && m!=k 
        return true
    else
        return false
    end
end

function velocity(state::WorldState)
    n_boids = state.n_boids
    for k in 1:n_boids
        scale = sqrt((state.accel[k][2])^2 + ((state.accel[k][1])^2)) / state.max_accel
        if sqrt((state.accel[k][2])^2 + ((state.accel[k][1])^2)) > state.max_accel
            state.accel[k] = (state.accel[k][1] / scale^2 , state.accel[k][2] / scale^2 )
        state.vel_vector[k] = (state.vel_vector[k][1] + state.accel[k][1], state.vel_vector[k][2] + state.accel[k][2])  
        end
    end
    return nothing
end

function cohesion(state::WorldState)
    position = state.boids
    n_boids = state.n_boids
    for k in 1:n_boids
        sum_x=0
        sum_y=0
        count = 0 
        for m in 1:n_boids
            if bros_detector(state, k,m) 
                sum_x = sum_x + position[m][1]
                sum_y = sum_y + position[m][2]
                count+=1 
            end
        end
    avg_x = sum_x / count 
    avg_y = sum_y / count 
    ax = avg_x - position[k][1]
    ay = avg_y - position[k][2]
    state.accel[k] = (state.accel[k][1] + ax * state.weight_coh, state.accel[k][2] + ay * state.weight_coh)
    end
    return nothing
end

function separation(state::WorldState)
    position = state.boids
    n_boids = state.n_boids
    for k in 1:n_boids 
        sum_x=0
        sum_y=0
        count = 0
        for m in 1:n_boids
            if bros_detector(state, k, m) 
                sum_x = sum_x + position[m][1]
                sum_y = sum_y + position[m][2]
                count+=1 
            end
        end
    avg_x = sum_x / count 
    avg_y = sum_y / count 
    ax = -(avg_x - position[k][1])
    ay = -(avg_y - position[k][2])
    state.accel[k] = (state.accel[k][1] + ax * state.weight_sep, state.accel[k][2] + ay * state.weight_sep)
    end
end
    

function alignment(state::WorldState)
    n_boids = state.n_boids
    for k in 1:n_boids 
        sum_x=0
        sum_y=0 
        count=0
        for m in 1:n_boids
            if bros_detector(state, k, m) 
                sum_x = sum_x + state.vel_vector[m][1]
                sum_y = sum_y + state.vel_vector[m][2]
                count+=1 
            end
        end
        avg_x = sum_x / count
        avg_y = sum_y / count
        ax = avg_x + state.vel_vector[k][1]
        ay = avg_y + state.vel_vector[k][2]
        state.accel[k] = (state.accel[k][1] + ax * state.weight_al, state.accel[k][2] + ay * state.weight_al)
    end
    return nothing
end


function maxspeed(state::WorldState)
    for k in 1:state.n_boids
    speed = sqrt(state.vel_vector[k][1]^2 + state.vel_vector[k][2]^2)
        if speed > state.max_vel
            scale = state.max_vel / speed
            state.vel_vector[k] = (state.vel_vector[k][1] * scale^2, state.vel_vector[k][2] * scale^2)
        end
    end 
    return nothing
end

function borders(state::WorldState)
    for k in 1:state.n_boids
        if state.boids[k][1] ≥ state.width
            state.vel_vector[k] = (-state.vel_vector[k][1], state.vel_vector[k][2])
            state.boids[k] = (state.width - 0.1, state.boids[k][2]) 
        elseif state.boids[k][1] ≤ 0
            state.vel_vector[k] = (-state.vel_vector[k][1], state.vel_vector[k][2])
            state.boids[k] = (0.1, state.boids[k][2]) 
        end

        if state.boids[k][2] ≥ state.height
            state.vel_vector[k] = (state.vel_vector[k][1], -state.vel_vector[k][2])
            state.boids[k] = (state.boids[k][1], state.height - 0.1) 
        elseif state.boids[k][2] ≤ 0
            state.vel_vector[k] = (state.vel_vector[k][1], -state.vel_vector[k][2])
            state.boids[k] = (state.boids[k][1],  0.1) 
        end
    end
    return nothing
end

function update!(state::WorldState)
    n_boids = state.n_boids
    for k in 1:n_boids
        borders(state) 
        state.boids[k] = state.boids[k] .+ state.vel_vector[k] 
        state.accel[k] = (0.0, 0.0)
    end
    cohesion(state)
    separation(state)
    alignment(state)
    velocity(state)
    maxspeed(state)
    return nothing 
end

function (@main)(ARGS)
    h = 10
    w = 10
    state = WorldState(h, w)
    anim = @animate for time = 1:100
        update!(state)
        boids = state.boids
        scatter(boids, xlim = (0, state.width), ylim = (0, state.height))
    end
    gif(anim, "boids.gif", fps =10)
end

end
#dew
using .Boids
Boids.main("")
