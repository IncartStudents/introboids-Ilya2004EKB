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
    max_vel::Float64
    closest_bro::Vector{Tuple{Float64, Float64}}
    function WorldState(n_boids, h, w)
        max_vel = 0.5
        # vel_module= [(rand() * max_vel ) for _ in 1:n_boids]
        vel_module= [(max_vel) for _ in 1:n_boids]
        boids = [(rand(0:w), rand(0:h)) for _ in 1:n_boids]
        vel_angle = [(rand() * 2π - π) for _ in 1:n_boids]
        vel_vector = [(vel_module[k] * cos(vel_angle[k]), vel_module[k] * sin(vel_angle[k])) for k in 1:n_boids]
        accel = [(0.0, 0.0) for _ in 1:n_boids]
        closest_bro = [(NaN, NaN) for _ in 1:n_boids]
        new(boids, vel_vector, vel_angle, vel_module, accel, h, w, max_vel, closest_bro)
    end
end

function bros_detector(state::WorldState, position, distance, k, m, sep_fact)
    x0=position[k][1]
    y0=position[k][2]
    x=position[m][1]
    y=position[m][2]
    if (abs(x-x0) ≤ distance) && (abs(y-y0) ≤ distance) && m!=k 
        if sqrt((x - x0)^2 + (y - y0)^2) ≤ sep_fact
            state.closest_bro[k]= (x0 - x, y0 - y)
        end
        return [x, y, state.vel_vector[m][1], state.vel_vector[m][2]]
    end
end

function velocity(state::WorldState, k, alignment, weight_al)
    state.vel_vector[k] = (state.vel_vector[k][1] + state.accel[k][1], state.vel_vector[k][2] + state.accel[k][2])
    # println(k, ' ',x, ' ', y, ' ',  state.vel_angle[k])
    return nothing
end

function cohesion(state::WorldState, position, grp_x, grp_y, k, max_accel)
    sum_x = sum(grp_x)
    sum_y = sum(grp_y)
    avg_x = sum_x / length(grp_x)
    avg_y = sum_y / length(grp_y)
    rx = avg_x - position[k][1]
    ry = avg_y - position[k][2]
    return (rx, ry)
end

function separation(state::WorldState, k, max_accel)
    if isnan(state.closest_bro[k][1]) == 0 
        x = state.closest_bro[k][1]
        y = state.closest_bro[k][2]
        return (x, y)
    end
end
    

function alignment(state::WorldState, k, grp_vel_x, grp_vel_y)
    sum_x = sum(grp_vel_x) 
    sum_y = sum(grp_vel_y) 
    avg_x = sum_x / length(grp_vel_x)
    avg_y = sum_y / length(grp_vel_x)
    adj_x = avg_x - state.vel_vector[k][1]
    adj_y = avg_y - state.vel_vector[k][2]
    return (adj_x, adj_y)
end

function compare(a,b,c,a_weight, b_weight,c_weight, max_accel)
    if isnothing(a) 
        a = (0.0 , 0.0)
    end
    if isnothing(b)
        b = (0.0 , 0.0)
    end
    if isnothing(c)
        c = (0.0 , 0.0)
    end
    x = (a[1] * a_weight + b[1] * b_weight + c[1] * c_weight) / 3
    y = (a[2] * a_weight + b[2] * b_weight + c[2] * c_weight) / 3
    accel_angle = atan(x/y)
    if (x < 0 && y < 0)
        accel_angle = -accel_angle - π/2
    elseif (x < 0 && y > 0) 
        accel_angle = -accel_angle + π/2
    end
    if sqrt(x^2 + y^2) > max_accel
        return (max_accel * cos(accel_angle), max_accel * sin(accel_angle))
    end
    return (x, y)
end

function maxspeed(state::WorldState, k)
    speed = sqrt(state.vel_vector[k][1]^2 + state.vel_vector[k][2]^2)
    if speed > state.max_vel
        scale = state.max_vel / speed
        state.vel_vector[k] = (state.vel_vector[k][1] * scale, state.vel_vector[k][2] * scale)
    end
    return nothing
end

function borders(state::WorldState, k)
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
    return nothing
end

function update!(state::WorldState, n_boids)
    state.closest_bro = [(NaN, NaN) for _ in 1:n_boids]
    for k in 1:n_boids
        borders(state, k) 
        state.boids[k] = state.boids[k] .+ state.vel_vector[k] 
    end

    distance = 6
    sep_fact = 2

    weight_coh = 0.4
    weight_sep = 2
    weight_al = 0.5

    max_accel = 0.5
    for k in 1:n_boids
        groups_x = Float64[]
        groups_y = Float64[]
        groups_vel_x =Float64[]
        groups_vel_y =Float64[]
        for m in 1:n_boids
            bros = bros_detector(state, state.boids, distance, k, m, sep_fact) 
            if bros !== nothing
                push!(groups_x, bros[1])
                push!(groups_y, bros[2])
                push!(groups_vel_x, bros[3])
                push!(groups_vel_y, bros[4])
            end
        end
        if isempty(groups_x) == 0
            coh = cohesion(state, state.boids, groups_x, groups_y, k, max_accel)
            sep = separation(state, k, max_accel)
            # println(coh, ' ', k)
            # println()
            al = alignment(state, k, groups_vel_x, groups_vel_y)
            acceleration = compare(coh, sep, al, weight_coh, weight_sep, weight_al, max_accel)
            state.accel[k] = acceleration 
            # state.accel[k] = (cohesion(state, state.boids, groups_x, groups_y, k, max_accel) .* weight_coh)
            velocity(state, k, al, weight_al)
            maxspeed(state, k)
        end
    end
    return nothing
end

function (@main)(ARGS)
    h = 50
    w = 50
    n_boids =  50
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
