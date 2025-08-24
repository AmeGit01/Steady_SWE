# FUNCTIONS for pseudo-transient method applied to 1D SWE with fixed bed

@views function build_bed!(z, x, iF, n_points)
    # build the bed of the channel, the idea is to have always positive z coordinates, and locate
    # the zero of the x coordinate at the most upstream point of the channel
    dz_tot = x[end] * iF                # total z coordinate at the end of the channel, 
                                        # it is a linear function of x
    for n in 1:n_points
        z[n] = -iF * x[n] + dz_tot      # z coordinate is a linear function of x, with slope iF
    end
    return nothing
end

@views function evaluate_depth(Q::Float64, B::Float64, Ks::Float64, iF::Float64)            # to be improoved
    # Find the normal depth from the formula above, it is an implicit formula, so adopt an iterative procedure:
    Qnew = 0.0;                                 # initialize discharge
    h = 1.0;                                    # initial guess for the depth [m]
    initial_res = Q - Q_formula(B, h, Ks, iF)   # save the initial residual
    max_iter = 10000                             # maximum number of iterations
    n = 0                                       # iteration counter
    tol = 1.0e-5; res = 0.1                     # tolerance and iniatalize the residual
    while abs(res) ≥ tol 
        Qnew = Q_formula(B, h, Ks, iF); res = Q - Qnew;     # evaluate discharghe with the new d
        increment = res/initial_res*0.01;                    # increment toward the solution
        if res > 0          # the depth is too small
            h += increment; # increase the depth
        else                # the depth is too high
            h -= increment; # decrease the depth
        end
        n += 1                                              # step
        # @printf("n = %d, d = %.5f \n", n, d)
        if n > max_iter                                     # control condition (avoid infinite cicle)
            error("Maximum number of iterations reached")
        end
    end
    c = Critical_d(Q, B)                                    # evaluate the critical depth
    return h
end 

Q_formula(B::Float64, h::Float64, Ks::Float64, iF::Float64) = B * h * Ks * Hyd_radius(B,h)^(2/3) * iF^(1/2);
Q_formula(B::Float64, h::Vector{Float64}, Ks::Float64, iF::Float64) = B .* h .* Ks .* Hyd_radius(B,h).^(2/3) .* iF^(1/2);
Hyd_radius(B::Float64, h::Float64) = (B*h)/(2*h + B)
Hyd_radius(B::Float64, h::Vector{Float64}) = (B.*h)./(2 .*h .+ B)
Critical_d(Q, B) = (Q^2/(9.81*B^2))^(1/3)