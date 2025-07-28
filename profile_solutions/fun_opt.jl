@views function evaluate_depths(goal::Float64, B::Vector{Float64}, Ks::Vector{Float64}, iF::Vector{Float64})
    N = length(iF)
    d_vec = Float64[]; c_vec = Float64[]
    max_iter = 100; tol = 1.e-5
    for n ∈ 1:N
        dL = 0.01; dR = 10; d = 0.5*(dL+dR) # meters
        QL = Q_formula(B[n], dL, Ks[n], iF[n]); δQL = QL - goal
        QR = Q_formula(B[n], dR, Ks[n], iF[n]); δQR = QR - goal
        if (δQL * δQR > 0) # check the initialization
            error("Increase the initial range of depths (bisection)")
        end
        Q  = Q_formula(B[n], d,  Ks[n], iF[n]); res = goal - Q;
        m = 0
        while abs(res) > tol
            m > max_iter && error("Maximum number of iterations reached (bisection)")
            m += 1
            if (res > 0)  # Q is too small, d must increase
                dL = d; QL = Q
            else
                dR = d; QR = Q
            end
            d = 0.5*(dL+dR); Q = Q_formula(B[n], d, Ks[n], iF[n])
            res = goal - Q
        end
        c = Critical_d(Q, B[n])
        
        append!(d_vec, d)
        append!(c_vec, c)
    end

    return d_vec, c_vec
end


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


function plot_1D(x, y, xlimits; lab = "undefined", xlab = "undefined", ylab = "undefined",
                 ttl = "undefined", cl = "undefined", grid = true) 
    # a flag for eventually limits of the plots..
    p1 = plot(x, y, label=lab, xlabel=xlab, ylabel=ylab,
    title=ttl, grid=true, color=cl, linewidth=2,
    xlim=xlimits)
    return p1
end


@views function analytical_E2d(d_critical, E, Q, B; style)
    # The ingredients are: 
    #  - Γ0 = E/Ec, 
    #  - Ec = critical energy
    #  - α  = arctan(sqrt(Γ0^3 - 1))
    # Equation to be solved:
    #  Γ0 - 2/3η - 1/3 η^-2 = 0, η = d/dc

    g  = PARAMETERS.gravit 
    Ec = Energy(Q, B, d_critical, g)                # Calculate critical energy (minimum)
    Γ0 = E/Ec                                       # Dimensionless energy
    @infiltrate false
    if Γ0 ≤ 1.0                                     # control E ≥ Ec
        η = 1.0                                     # if not, d = d_criical
    else
        α  = atan(sqrt(Γ0^3 - 1))                   # coefficient
        if style == "supercritical"
            η = 0.5*Γ0*(1+2*cos((π+2*α)/3)) # root 3, dimensionless supercritical depth
        elseif style == "subcritical"
            η = 0.5*Γ0*(1+2*cos((π-2*α)/3)) # root 2, dimensionless subcritical depth
        else
            error("Style must be either 'subcritical' or 'supercritical'")
        end
    end
    # roots 2 and 3 are given from: 
    # [ ] A. Viviani and V. Caleffi, Depth-energy and depth-force relationship 
    # in open channel flows: Analytical findings, Advances in Water Resouces 31, Elsevier, 2008
    d = η * d_critical                              # dimensional depth
    return d
end

Spinta(Q, B, d) = 500*9.81*B*d^2 + 1000*Q^2/(B*d)
Energy(Q, B, d, g) = d + Q^2/(2*g*(B*d)^2)
Hyd_radius(B, d) = (B*d)/(2*d + B)
dEdx(iF, Q, B, d, Ks) = iF - ( Q^2 / (Ks^2 * B^2 * d^2) * Hyd_radius(B,d)^(-4/3) )
Q_formula(B, d, Ks, iF) = B * d * Ks * Hyd_radius(B,d)^(2/3) * iF^(1/2);
Critical_d(Q, B) = (Q^2/(9.81*B^2))^(1/3)