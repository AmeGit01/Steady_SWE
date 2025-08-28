# Main code fo pseudo-transient method applied to 1D SWE with fixed bed
using IJulia, Plots, Printf, LinearAlgebra, BenchmarkTools, Infiltrator
include("functions_PT_1D.jl")


# So far the idea is to try the method for a single reach case, 
# in which slope and width are constant, and also the inflow discharge
# The initial condition will be a fixed water depth different than the 
# elevation at the steady solution, to see whether the method converges or not
# The boundary condition at left will be the discharge, and at right the water depth

const PARAMETERS = (
    gravit = 9.81,     # m/s^2
    ρ      = 1000.0   # kg/m^3
)

@views function Pseudotransient_1D()
    # Physical parameters
    Q = 20.0                # discharge [m^3/s]
    B = 10.0                # channel width [m]
    Ks = 30.0               # Strickler coefficient [m^(1/3)/s]
    iF = 0.001              # channel slope [-]
    L = 1000.0              # channel length [m]
    n_points = 10           # number of points
    g = PARAMETERS.gravit
    ρ = PARAMETERS.ρ
    γ = g * ρ               # specific weight of water [N/m^3]

    # Numerical parameters
    dx = L/(n_points)                   # space step [m]
    xb = collect(dx/2:dx:L-dx/2)        # coordinate of the barycenters  (n_points)
    xi = collect(0:dx:L)                # coordinate of the interfaces (n_points+1)
    ρ_hat_bar = 1.0                     # parameter defining the pseudo-transient method
    tol = 1.e-8                         # tolerance for the convergence
    max_iter = 1000                     # maximum number of iterations

    # Initialize arrays
    b     = zeros(n_points)                # bed elevation [m] BARYCENTERS
    h     = ones(n_points)                 # water depth [m]  BARYCENTERS (INITIAL VALUES ALREADY INITIALIZED)
    q     = zeros(n_points+1)
    u     = zeros(n_points+1)
    iF_v  = ones(n_points+1) * iF
    hm    = zeros(n_points+1)
    hgost = zeros(n_points+1)
    dudx  = zeros(n_points)

    # Build the bed
    build_bed!(b, xb, iF, n_points)

    # Evaluate the normal depth for the selected discharge
    hL = evaluate_depth(Q, B, Ks, iF)

    # set the initial conditions
    h[1] = hL; h_old = copy(h)
    hm[1] = h_old[1]; hm[end] = h_old[end]
    hm[2:end-1] .= 0.5 .* (h_old[1:end-1] .+ h_old[2:end])

    q .= Q_formula(B, hm, Ks, iF); q_old = copy(q)
    # q .= Q;    q_old = copy(q)
    # u .= Q./h

    iter = 1; err = 2tol
    while err ≥ tol && iter ≤ max_iter
        h[2:end] .= h[2:end] .- ρ / ρ_hat_bar / dx .* diff(q[2:end])
        @infiltrate false
        hm[1] = h[1]; hm[end] = h[end]
        hm[2:end-1] .= 0.5 .* (h[1:end-1] .+ h[2:end])
        u .= q ./ hm
        dudx .= diff(u) ./ dx
        hgost[1:end-1] .= h_old; hgost[end] = h_old[end] 
        q[2:end] .= q[2:end] .- 1/ρ_hat_bar .*( ρ .* q[2:end] .* dudx .-
                    γ .* hm[2:end] .* ( iF_v[2:end] .- diff(hgost) .- 
                    u[2:end].^2 ./(Ks^2 .* hm[2:end].^4 .* hm[2:end].^(-3)) ) )
        @infiltrate false
        
        any(h==NaN) && error("h is NaN")
        any(q==NaN) && error("q is NaN")
        # Convergence check
        err = max(maximum(abs.(h .- h_old)), maximum(abs.(q .- q_old)))
        @printf("iter = %d, err = %.5e \n", iter, err)

        # Update the solutions
        h_old .= h; q_old .= q
        iter += 1
    end
end

Pseudotransient_1D()


#=  DO NOT WORK, BUT COULD BE A STARTING POINT
    iter = 1; err = 2tol
    while err ≥ tol && iter ≤ max_iter
        h[2:end] .= h_old[2:end] .- ρ / ρ_hat_bar / dx .* diff(q_old[2:end])
        @infiltrate false
        hm[1] = h_old[1]; hm[end] = h_old[end]
        hm[2:end-1] .= 0.5 .* (h_old[1:end-1] .+ h_old[2:end])
        dudx .= diff(q_old./hm) ./ dx
        hgost[1:end-1] .= h_old; hgost[end] = h_old[end] 
        u .= q_old ./ hm
        q[2:end] .= q_old[2:end] .- 1/ρ_hat_bar .*( ρ .* q_old[2:end] .* dudx .-
                    γ .* hm[2:end] .* ( iF_v[2:end] .- diff(hgost) .- 
                    u[2:end].^2 ./(Ks^2 .* hm[2:end].^4 .* hm[2:end].^(-3)) ) )
        @infiltrate false
        
        any(h==NaN) && error("h is NaN")
        any(q==NaN) && error("q is NaN")
        # Convergence check
        err = max(maximum(abs.(h .- h_old)), maximum(abs.(q .- q_old)))
        @printf("iter = %d, err = %.5e \n", iter, err)

        # Update the solutions
        h_old .= h; q_old .= q
        iter += 1
    end
=#