@views function evaluate_depths(Q::Float64, B::Float64, Ks::Float64, iF::Float64)
    # Find the normal depth from the formula above, it is an implicit formula, so adopt an iterative procedure:
    Qnew = 0.0;                                 # initialize discharge
    d = 1.0;                                    # initial guess for the depth [m]
    initial_res = Q - Q_formula(B, d, Ks, iF)   # save the initial residual
    max_iter = 10000                             # maximum number of iterations
    n = 0                                       # iteration counter
    tol = 1.0e-5; res = 0.1                     # tolerance and iniatalize the residual
    while abs(res) ≥ tol 
        Qnew = Q_formula(B, d, Ks, iF); res = Q - Qnew;     # evaluate discharghe with the new d
        increment = res/initial_res*0.01;                    # increment toward the solution
        if res > 0          # the depth is too small
            d += increment; # increase the depth
        else                # the depth is too high
            d -= increment; # decrease the depth
        end
        n += 1                                              # step
        # @printf("n = %d, d = %.5f \n", n, d)
        if n > max_iter                                     # control condition (avoid infinite cicle)
            error("Maximum number of iterations reached")
        end
    end
    c = Critical_d(Q, B)                                    # evaluate the critical depth
    return d, c
end 

# Vectorial form
@views function evaluate_depths(Q::Float64, B::Vector{Float64}, Ks::Vector{Float64}, iF::Vector{Float64})
    N = length(IF)
    d_vec = Float64[]; c_vec = Float64[]
    for n ∈ 1:N
        # Find the normal depth from the formula above, it is an implicit formula, so adopt an iterative procedure:
        Qnew = 0.0;                                 # initialize discharge
        d = 1.0;                                    # initial guess for the depth [m]
        initial_res = Q - Q_formula(B[n], d, Ks[n], iF[n])   # save the initial residual
        max_iter = 10000                             # maximum number of iterations
        m = 0                                       # iteration counter
        tol = 1.0e-5; res = 0.1                     # tolerance and iniatalize the residual
        while abs(res) ≥ tol 
            Qnew = Q_formula(B[n], d, Ks[n], iF[n]); res = Q - Qnew;     # evaluate discharghe with the new d
            @infiltrate false
            increment = res/initial_res*0.01;                    # increment toward the solution
            if res > 0          # the depth is too small
                d += increment; # increase the depth
            else                # the depth is too high
                d -= increment; # decrease the depth
            end
            m += 1                                              # step
            if m > max_iter                                     # control condition (avoid infinite cicle)
                error("Maximum number of iterations reached")
            end
            # @printf("iter = %d \n", m)
        end
        c = Critical_d(Q, B[n])                                    # evaluate the critical depth

        append!(d_vec, d)
        append!(c_vec, c)
    end

    return d_vec, c_vec
end 

# scalar form with bisection
@views function evaluate_depths_b(goal::Float64, B::Float64, Ks::Float64, iF::Float64)
    # same as before but using bisection method (convergence guaranteed)
    dL = 0.01; dR = 10; d = 0.5*(dL+dR) # meters
    QL = Q_formula(B, dL, Ks, iF); δQL = QL - goal
    QR = Q_formula(B, dR, Ks, iF); δQR = QR - goal
    if (δQL * δQR > 0) # check the initialization
        error("Increase the initial range of depths (bisection)")
    end
    Q  = Q_formula(B, d,  Ks, iF); res = goal - Q;
    max_iter = 100; tol = 1.e-5; n = 0
    while abs(res) > tol
        n > max_iter && error("Maximum number of iterations reached (bisection)")
        n += 1
        if (res > 0)  # Q is too small, d must increase
            dL = d; QL = Q
        else
            dR = d; QR = Q
        end
        d = 0.5*(dL+dR); Q = Q_formula(B, d, Ks, iF)
        res = goal - Q
    end
    c = Critical_d(Q, B) 
    return d, c
end

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
                n += 1
                if n > max_iter
                    return d = 0.0
                end
            end
        elseif (style=="supercritical")
            while abs(res) ≥ 1.e-6
                @infiltrate false
                if res < 0
                    dL = d; fL = f;
                else
                    dR = d; fR = f;
                end
                d = (dL + dR) / 2.0; f = Energy(Q, B, d, g)
                res = goal - f
                n += 1
                if n > max_iter
                    return d = 0.0
                end
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

@views function analitical_E2d(d_critical, E, Q, B, g, n; style)
    # The ingredients are: 
    #  - Γ0 = E/Ec, 
    #  - Ec = critical energy
    #  - α  = arctan(sqrt(Γ0^3 - 1))
    # Equation to be solved:
    #  Γ0 - 2/3η - 1/3 η^-2 = 0, η = d/dc

    Ec = Energy(Q, B, d_critical, g)
    Γ0 = E/Ec
    if (n==967)
        @infiltrate false
    end
    
    if Γ0 ≤ 1.0
        η = 0.0
    else
        α  = atan(sqrt(Γ0^3 - 1))
        if style == "supercritical"
            η = 0.5*Γ0*(1+2*cos((π+2*α)/3)) # root 3
        elseif style == "subcritical"
            η = 0.5*Γ0*(1+2*cos((π-2*α)/3)) # root 2
        else
            error("Style must be either 'subcritical' or 'supercritical'")
        end
    end
    d = η * d_critical
    return d
end

@views function main_E_solver(Q, N, B_v, KS, IF)

    # I have to make the code working also for changing bed elevation, so the 
    # "position" of the reaches can change, i.e. the position of the passage 
    # through the critical state change avery time, so the code must be able to 
    # recognize itself where and how many they are. 

    # The only thing that remain constant for all the simulation is the 
    # total number of points (N) and the X coordinate (and dX vector), so I have to 
    # refer to them for doing all the work.

    # --------------------- BOUNDARY CONDITIONS ------------------------
    g = PARAMETERS.gravit

    @infiltrate false
    # vectorial version
    d_uniform = d_critical = zeros(N)
    d_uniform, d_critical = evaluate_depths(Q, B_v, KS, IF)

    @printf("Critical depth = %.3f m \n", d_critical[Int(floor(N/2))])

    # check wether the reaches are subcritical or supercritical
    subcritical = zeros(N)              # vector describing wether the reach is subcritical or supercritical
    n_changes = 0                       # counter for the changes in style
    ch_loc = Int64[]                    # location of the changes, needed for the BC
    for n ∈ 1:N
        if d_uniform[n] > d_critical[n]
            subcritical[n] = 1
        else 
            subcritical[n] = 0
        end
        # count and print the number of changes in style
        if n>1 && subcritical[n] ≠ subcritical[n-1]
            n_changes += 1
            append!(ch_loc, n)
        end 
    end
    n_reaches = n_changes + 1
    @infiltrate false
    @printf("Number of style changes is %d so there are %d reaches. \n", n_changes, n_reaches)


    # evaluate the BC for each reach (considered separately), so that to impose them internally where requested
    EL = zeros(n_reaches); dL = zeros(n_reaches) # left boundary condition for the energy
    ER = zeros(n_reaches); dR = zeros(n_reaches) # right boundary condition for the energy

    EL[1]         = Energy(Q, B_v[1],         d_uniform[1],         g); 
    dL[1]         = d_uniform[1]  # water elevation at the left boundary [m]

    ER[n_reaches] = Energy(Q, B_v[end], d_uniform[end], g); 
    dR[n_reaches] = d_uniform[end] # water elevation at the right boundary [m]

    m = 2
    for n_rea in ch_loc
        discriminant = subcritical[n_rea] - subcritical[n_rea-1] # check if the reach is subcritical or supercritical
        if discriminant == 0 || discriminant == +1
            EL[m]   = Energy(Q, B_v[n_rea],   d_uniform[n_rea],   g); dL[m]   = d_uniform[n_rea] # if the reach is subcritical, use the uniform depth
            ER[m-1] = Energy(Q, B_v[n_rea-1], d_uniform[n_rea-1], g); dR[m-1] = d_uniform[n_rea]
        elseif discriminant == -1
            EL[m]   = Energy(Q, B_v[n_rea],   d_critical[n_rea], g);   dL[m]   = d_critical[n_rea] # if the reach is supercritical, use the critical depth
            ER[m-1] = Energy(Q, B_v[n_rea-1], d_critical[n_rea-1], g); dR[m-1] = d_critical[n_rea-1] # if the previous reach is supercritical, use the critical depth
        else
            error("Discriminant must be either 0, -1, or 1")
        end
        m += 1
    end

    # Print the BC for checking purposes
    for n_rea in 1:n_reaches
        @printf("Reach %d: Left E BC = %.4f m, Right E BC = %.4f m \n", n_rea, EL[n_rea], ER[n_rea])
        @printf("Reach %d: Left d BC = %.4f m, Right d BC = %.4f m \n", n_rea, dL[n_rea], dR[n_rea])
    end 

    # ------------------------- INTEGRATION ----------------------------- 
    # NB: here we have to solve for E and immediately evaluate the depths from the energy

    # Solve the energy equation from left to right (downward)
    # here we have to solve the supercritical reaches, starting from the left boundary and
    # going on until reaching the critical state
    e_dw   = zeros(N)
    d_fast = zeros(N) # depth from the downward energy
    ch_loc_dw = copy(ch_loc)
    insert!(ch_loc_dw, 1, 1)
    m = 0
    for n_rea ∈ ch_loc_dw
        m += 1 
        # check the reach is supercritical
        if subcritical[n_rea] == 1
            @printf("Reach %d is subcritical, skipping the supercritical computation \n", n_rea)
            continue # skip the supercritical computation for subcritical reaches
        elseif n_rea>1 && (subcritical[n_rea-1] == 0 && subcritical[n_rea] == 0)
            continue # skip the computation if the reach is supercritical and the previous one was supercritical too
        else
            @printf("Supercritical computation for reach %d \n", m)
        end

        # set the beginning of the reach
        n = n_rea + 1
        e_dw[n] = EL[m]
        d_fast[n] = dL[m]
        @printf("Imposed BC: depth = %.4f m,  energy = %.4f m \n", d_fast[n], e_dw[n])

        while n < N # evaluate energy for downstream points
            n += 1
            # Forward Euler
            e_dw[n] = e_dw[n-1] + dX[n-1] * dEdx(IF[n-1], Q, B_v[n-1], d_fast[n-1], KS[n-1])
            # evaluate uniform depth (supercritical)
            d_fast[n] = analytical_E2d(d_critical[n], e_dw[n], Q, B_v[n], style="supercritical")
            # when the criticality is reached, than stop the evaluation:
            @infiltrate false
            d_fast[n] == d_critical[n] && break

            # @printf("n = %d, e_dw[n] = %f, d_fast[n] = %.4f, X[n] = %.2f \n", n, e_dw[n], d_fast[n], X[n])
        end
    end

    y_limits = (min(minimum(Z), minimum(e_dw)) - 0.5, max(maximum(Z), maximum(Z+d_fast)) + 0.5) # y limits for the plots
    p1 = plot(X, Z, label="Bed", xlabel="x [m]", ylabel="z [m]", title="Bed of the channel", legend=:topright,
        grid=true, xlims=(X[begin]-10, X[end]+10), ylims=y_limits, color="black", linewidth=2,
        size=(1400, 700))
    p1 = plot!(X, Z + d_fast, label="Supercritical depths", xlabel="x [m]", ylabel="z [m]",
        title="Water elevation in the channel", grid=true, color="red", linewidth=2)
    p1 = plot!(X, Z + d_critical, label="Critical depth", xlabel="x [m]", ylabel="z [m]",
        title="Water elevation in the channel", grid=true, color="gray", linewidth=2, linestyle=:dash)
    display(plot(p1, size=(1400, 500)))


    # Solve the energy equation from right to left (upward)
    # here we have to solve the subcritical reaches, starting from the right boundary and
    # going on until reaching the critical state
    e_uw   = zeros(N) 
    d_slow = zeros(N) 
    ch_loc_uw = copy(ch_loc)
    append!(ch_loc_uw, N+1)
    ch_loc_uw .-= 1
    reverse!(ch_loc_uw)
    m = length(ch_loc_uw)+1
    for n_rea ∈ ch_loc_uw
        m -= 1
        # check the reach is subcritical
        if subcritical[n_rea] == 0
            @printf("Reach %d is supercritical, skipping the subcritical computation \n", n_rea)
            continue # skip the supercritical computation for subcritical reaches
        elseif n_rea < N && (subcritical[n_rea+1] == 1 && subcritical[n_rea] == 1)
            continue
        else
            @printf("Subcritical computation for reach %d \n", m)
        end

        # set the end of the reach
        n = n_rea
        e_uw[n] = ER[m]
        d_slow[n] = dR[m] # initial guess for the depth
        @printf("Imposed BC: depth = %.4f m,  energy = %.4f m \n", d_slow[n], e_uw[n])

        while n > 1
            n -= 1
            # Forward Euler (from dowstream to upstream)
            e_uw[n] = e_uw[n+1] - dX[n] * dEdx(IF[n+1], Q, B_v[n+1], d_slow[n+1], KS[n+1])
            # evaluate uniform depth (subcritical)
            d_slow[n] = analytical_E2d(d_critical[n], e_uw[n], Q, B_v[n], style="subcritical")
            # when the criticality is reached, than stop the evaluation:
            d_slow[n] == d_critical[n] && break

            # @printf("n = %d, e_uw[n] = %f, d_slow[n] = %.4f, X[n] = %.2f \n", n, e_uw[n], d_slow[n], X[n])
        end
    end

    # ------------------------- SPINTA EVALUATION -----------------------

    # evaluate the Spinta for both the supercritical and subcritical depths
    Spinta_fast = zeros(N)                   # Spinta for the supercritical depths  
    Spinta_slow = zeros(N)                   # Spinta for the subcritical depths
    d           = zeros(N)                   # final depth to be plotted
    for n = 1:N
        Spinta_fast[n] = d_fast[n]==0 ? 0 : Spinta(Q, B_v[n], d_fast[n])    # Spinta for the supercritical depths
        Spinta_slow[n] = d_slow[n]==0 ? 0 : Spinta(Q, B_v[n], d_slow[n])    # Spinta for the subcritical depths
        # use the supercritical depth if Spinta is greater, otherwise use the subcritical depth
        Spinta_fast[n] > Spinta_slow[n] ? d[n] = d_fast[n] : d[n] = d_slow[n] 
    end
    #println("Spinta for the supercritical depths: ", Spinta_fast[1], " N/m, with depth: ", d_fast[1], " m")
    #println("Spinta for the subcritical depths: ", Spinta_slow[1], " N/m, with depth: ", d_slow[1], " m")

    # --------------------------- PLOTS -----------------------------

    y_limits = (min(minimum(Z), minimum(e_dw)) - 0.5, max(maximum(Z), maximum(Z+d)) + 0.5) # y limits for the plots
    p1 = plot(X, Z, label="Bed", xlabel="x [m]", ylabel="z [m]", title="Bed of the channel", legend=:topright,
        grid=true, xlims=(X[begin]-10, X[end]+10), ylims=y_limits, color="black", linewidth=2,
        size=(1400, 700))
    p1 = plot!(X, Z + d_fast, label="Supercritical depths", xlabel="x [m]", ylabel="z [m]",
        title="Water elevation in the channel", grid=true, color="red", linewidth=2)
    p1 = plot!(X, Z + d_slow, label="Subcritical depths", xlabel="x [m]", ylabel="z [m]",
        title="Water elevation in the channel", grid=true, color="blue", linewidth=2)
    p1 = plot!(X, Z + d, label="Water depth", xlabel="x [m]", ylabel="z [m]",
        title="Water elevation in the channel", grid=true, color="green", linewidth=2)
    p1 = plot!(X, Z + d_critical, label="Critical depth", xlabel="x [m]", ylabel="z [m]",
        title="Water elevation in the channel", grid=true, color="gray", linewidth=2, linestyle=:dash)
    display(plot(p1, size=(1400, 500)))

    #=
    p2 = plot(X, Z, label="Bed", xlabel="x [m]", ylabel="z [m]", title="Bed of the channel", legend=:topright,
        grid=true, xlims=(950, 3050), ylims=(Z[3000]-1, Z[1000]+d_slow[100]+2), color="black", linewidth=2,
        size=(1400, 700))
    p2 = plot!(X, Z + d_fast, label="Supercritical depths", xlabel="x [m]", ylabel="z [m]",
        title="Water elevation in the channel", grid=true, color="red", linewidth=2)
    p2 = plot!(X, Z + d_slow, label="Subcritical depths", xlabel="x [m]", ylabel="z [m]",
        title="Water elevation in the channel", grid=true, color="blue", linewidth=2)
    p2 = plot!(X, Z + d, label="Water depth", xlabel="x [m]", ylabel="z [m]",
        title="Water elevation in the channel", grid=true, color="green", linewidth=2)
    p2 = plot!(X, Z + d_critical, label="Critical depth", xlabel="x [m]", ylabel="z [m]",
        title="Water elevation in the channel", grid=true, color="gray", linewidth=2, linestyle=:dash)
    =#
    # display(plot(p1, p2, layout=(2,1), size=(1400, 800)))
    
    return nothing
end #main_E_solver

Spinta(Q, B, d) = 500*9.81*B*d^2 + 1000*Q^2/(B*d)
Energy(Q, B, d, g) = d + Q^2/(2*g*(B*d)^2)
Hyd_radius(B, d) = (B*d)/(2*d + B)
dEdx(iF, Q, B, d, Ks) = iF - ( Q^2 / (Ks^2 * B^2 * d^2) * Hyd_radius(B,d)^(-4/3) )
Q_formula(B, d, Ks, iF) = B * d * Ks * Hyd_radius(B,d)^(2/3) * iF^(1/2);
Critical_d(Q, B) = (Q^2/(9.81*B^2))^(1/3)