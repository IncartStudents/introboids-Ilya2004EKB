module Boids
 import Pkg; Pkg.add("Plots")
using Plots

mutable struct WorldState
    boids::Vector{Tuple{Float64, Float64}}
    vel_vector::Vector{Tuple{Float64, Float64}}
    height::Float64
    width::Float64
    function WorldState(n_boids, h, w)
        boids = [(rand(0:w), rand(0:h)) for _ in 1:n_boids]
        for k in 1:n_boids
            vel_module =rand(0:1) * 0.1
            vel_angle = rand(-π:π)
            vel_vector =  [vel_module * cos(vel_angle), vel_module * sin(vel_angle)]
        end
        new(boids, h, w, vel_vector)
    end
end

function update!(state::WorldState, n_boids)
    for k in 1:n_boids
        state.boids[k] = broadcast(+, state.boids[k], state.vel_vector[k])
    end
    # TODO: 
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
