@views function evaluate_depths(Q, B, L, Ks, iF)
    # Find the normal depth from the formula above, it is an implicit formula, so adopt an iterative procedure:
    Qnew = 0.0; tol = 1.0e-4; res = 0.1; d = 1.0; # initial guess for the depth [m]
    initial_res = Q - B * d * Ks * ((B * d)/(2*d + B))^(2/3) * iF^(1/2);
    max_iter = 1000 # maximum number of iterations
    n = 0
    while abs(res) >= tol 
        Qnew = B * d * Ks * ((B * d)/(2*d + B))^(2/3) * iF^(1/2)
        res = Q - Qnew;
        increment = res/initial_res*0.1;
        if res > 0 # so the depth is too small
            d += increment; # increase the depth
        else
            d -= increment; # decrease the depth
        end
        n += 1 
        if n > max_iter
            error("Maximum number of iterations reached")
        end
    end
    c = (Q^2/(9.81*B^2))^(1/3)
    return d, c
end #function

@views function build_bed!(z, x, dx, iF, n_points)
    # build the bed of the channel, the idea is to have always positive z coordinates, and locate
    # the zero of the x coordinate at the most upstream point of the channel
    dz_tot = x[end] * iF            # total z coordinate at the end of the channel, it is a linear function of x
    for n in 1:n_points
        z[n] = -iF * x[n] + dz_tot # z coordinate is a linear function of x, with slope iF
    end
    # return z
end

@views function E2d(d_c, Q, B, e; style)
    g = 9.81 # gravity acceleration [m/s^2]
    if style == "subcritical" 
        d = bisection_function(d_c, d_c + 10.0, Q, B, g, e, style)
    elseif style == "supercritical"
        d = bisection_function(0.001, d_c, Q, B, g, e, style)
    else
        error("Style must be either 'subcritical' or 'supercritical'")
    end
    return d
end

@views function bisection_function(dL, dR, Q, B, g, goal, style)
    d = (dL + dR) / 2.0; f  = Energy(Q, B, d, g)
    fL = Energy(Q, B, dL, g); ffL = fL - f;
    fR = Energy(Q, B, dR, g); ffR = fR - f;
    n = 0; max_iter = 100
    if (ffL) * (ffR) > 0
        error("Function values at the endpoints must have opposite signs")
    else
        res = goal - f
        if (style=="subcritical")
            while abs(res) ≥ 1.e-6
                if res > 0
                    dL = d; fL = f;
                else
                    dR = d; fR = f;
                end
                d = (dL + dR) / 2.0; f = Energy(Q, B, d, g)
                res = goal - f
            end
        elseif (style=="supercritical")
            while abs(res) ≥ 1.e-6
                if res < 0
                    dL = d; fL = f;
                else
                    dR = d; fR = f;
                end
                d = (dL + dR) / 2.0; f = Energy(Q, B, d, g)
                res = goal - f
            end
        end
    end
    return d
end

@views function E2d_single_reach(d_c, Q, B, e; style)
    g = 9.81 # gravity acceleration [m/s^2]
    n = 0; max_iter = 1000
    if style == "subcritical" 
        d = d_c + 1.0
        obj = d + Q^2 / (2 * g * (B * d)^2)
        initial_res = obj - e; res = initial_res
        while abs(res) >=1.e-5
            obj = d + Q^2 / (2 * g * (B * d)^2)
            res = obj - e
            increment = res/initial_res*0.1
            if res > 0 # so the energy is too much
                d -= increment; # decrease the depth
            else
                d += increment; # increase the depth
            end
            n += 1 
            if n > max_iter
                error("Maximum number of iterations reached")
            end
        end
    elseif style == "supercritical"
        d = d_c - 1.0
        obj = d + Q^2 / (2 * g * (B * d)^2)
        initial_res = obj - e; res = initial_res
        while abs(obj - e) >=1.e-5
            obj = d + Q^2 / (2 * g * (B * d)^2)
            res = obj - e
            increment = res/initial_res*0.1
            if res > 0 # so the energy is too much
                d += increment; # increase the depth
            else
                d -= increment; # decrease the depth
            end
            n += 1 
            if n > max_iter
                error("Maximum number of iterations reached")
            end
        end
    else
        error("Style must be either 'subcritical' or 'supercritical'")
    end
    return d
end

function plot_1D(x, y, xlimits; lab = "undefined", xlab = "undefined", ylab = "undefined",
                 ttl = "undefined", cl = "undefined", grid = true) 
    # a flag for eventually limits of the plots..
    p1 = plot(x, y, label=lab, xlabel=xlab, ylabel=ylab,
    title=ttl, grid=true, color=cl, linewidth=2,
    xlim=xlimits)
end

Spinta(Q, B, d) = 500*9.81*B*d^2 + 1000*Q^2/(B*d)
Energy(Q, B, d, g) = d + Q^2/(2*g*(B*d)^2)