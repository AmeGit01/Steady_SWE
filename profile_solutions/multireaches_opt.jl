using IJulia, Plots, Printf, Infiltrator, Revise, Test
include("fun_opt.jl")
GC.gc()

# --------------------- PARAMETERS DEFINITION ---------------------
const PARAMETERS = (
    gravit = 9.81,
    useless = 1.0
)


function bed_construction!(B_v, KS, IF, X, Z, L, B, Ks, iF, dx)
    n_reaches = length(L)
    # -------------------------- GEOMETRY -----------------------------

    n_points = Array{Int64}(undef, n_reaches)                                       # number of points in the channel
    for n_rea in 1:n_reaches
        # Check the lengths of the channel are multiple of dx
        L[n_rea] % dx ≠ 0 && error("The length of the channel must be a multiple of dx")
        n_points[n_rea] = Int(L[n_rea]/dx) + 1
        @printf("Number of points along reach %d: %d \n", n_rea, n_points[n_rea])
    end

    n_points_total = sum(n_points)                                                  # total number of points in the channel
    @printf("Total number of points in the channel: %d \n", n_points_total)

    dX = dx .* ones(n_points_total-1)                       # vector of distances between the points
    for n_rea in 1:n_reaches-1
        location = sum(n_points[1:n_rea]) + 1               # location of the first point of the next reach
        dX[location] = 0.0                                  # set the last distance to zero, to avoid problems in the integration
    end


    for n_rea in 1:n_reaches
        append!(B_v, B[n_rea].*ones(n_points[n_rea]))      # vector of widths
        append!(KS, Ks[n_rea].*ones(n_points[n_rea]))      # vector of Strickler coefficients
        append!(IF, iF[n_rea].*ones(n_points[n_rea]))      # vector of slopes
    end

    println("size B_v =  ", size(B_v), ", size KS = ", size(KS), ", size IF = ", size(IF), ", size dX = ", size(dX))

    # Check the dX vector is corrected (must be present some zeros)
    pdx = scatter(dX, label="dX", xlabel="Point index", ylabel="dX [m]",
        title="Distance between the points", grid=true, color="black",
        xlims=(1, n_points_total), size=(1000, 400)); # display(pdx)    
    # savefig("pdx_plot.png")

    # auxiliary vectors
    scale = zeros(n_reaches, 2)
    for n_rea in n_reaches-1:-1:1
        scale[n_rea, 2] = scale[n_rea+1, 2] + L[n_rea+1] * iF[n_rea+1]        # scale vector for the z coordinates
    end
    for n_rea in 2:n_reaches
        scale[n_rea, 1] = scale[n_rea-1, 1] + L[n_rea-1]                      # scale vector for the x coordinates
    end

    println("Scale vector for X: ", scale[:, 1])
    println("Scale vector for Z: ", scale[:, 2])

    for n_rea in 1:n_reaches
        x = LinRange(0, L[n_rea], n_points[n_rea])                  # x coordinates of the points
        z = zeros(n_points[n_rea])                                  # initialize z coordinates of the points
        build_bed!(z, x, iF[n_rea], n_points[n_rea])            # build the bed of the channel
        append!(X, x .+ scale[n_rea, 1])
        append!(Z, z .+ scale[n_rea, 2])
    end

    println("first and last x coordinates: ", X[1], " m, ", X[end], " m")
    println("first and last z coordinates: ", Z[1], " m, ", Z[end], " m")

    pbed = plot(X, Z, label="Bed", xlabel="x [m]", ylabel="z [m]", title="Bed of the channel",
        grid=true, xlims=(X[begin]-10, X[end]+10), ylims=(minimum(Z)-0.05, maximum(Z)+0.05),
        color="black", linewidth=2, size=(1000, 400)); #display(pbed)


    px   = plot(X, label="X coordinates", xlabel="Point index", ylabel="X [m]",
        title="X coordinates of the points", grid=true, color="blue",
        xlims=(1, n_points_total), size=(1000, 400))
    # display(plot(pbed, px, layout=(2,1), size=(1000, 800), legend=:topright)) 
    return dX, n_points_total
end # bed_construction


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


# ---------------------- PROBLEM DEFINITION -----------------------
n_reaches = 3               # number of reaches in the channel
L  = [500.0, 1000.0, 8000.0, 500.0]       # length of the channel [m]
Q  = 60.0                  # discharge [m^3/s] Q = Ω Ks R^2/3 iF^1/2
B  = [50.0, 50.0, 50.0, 50.0]           # width of the channel [m]
Ks = [40.0, 40.0, 50.0, 40.0]           # Strickler coefficient [m^(1/3)/s] Ks = 1/n
iF = [0.001, 0.01, 0.0001, 0.001]          # slope of the channel
dx = 1.0                   # distance between the points of the channel [m]
g  = PARAMETERS.gravit             # gravity acceleration [m/s^2]

B_v = Float64[]                                        # vector of widths
KS  = Float64[]                                        # vector of Strickler coefficients
IF  = Float64[]                                        # vector of slopes
X   = Float64[]
Z   = Float64[]
dX  = Float64[]

# dX, N = bed_construction!(B_v, KS, IF, X, Z, L, B, Ks, iF, dx)

# main_E_solver(Q, N, B_v, KS, IF)

# test the bisection function:
Q  = [20.0, 50.0, 100.0, 200.0]
B  = [10.0, 10.0, 50.0, 50.0]
Ks = [20.0, 20.0, 30.0, 40.0]
iF = [0.01, 0.001, 0.01, 0.001]

# d = c = zeros(4)
d_solutions = Float64[]
c_solutions = Float64[]

for i ∈ 1:4
    local d, c = evaluate_depths(Q[i], B[i], Ks[i], iF[i])
    # @printf("i = %d, uniform = %.5f, critical = %.5f \n", i, d, c)
    append!(d_solutions, d)
    append!(c_solutions, c)
end


@testset "evaluate_depths_new" begin
    @test evaluate_depths_new(Q[1], B[1], Ks[1], iF[1]) ≈ (d_solutions[1], c_solutions[1])
    @test evaluate_depths_new(Q[2], B[2], Ks[2], iF[2]) ≈ (d_solutions[2], c_solutions[2])
    @test evaluate_depths_new(Q[3], B[3], Ks[3], iF[3]) ≈ (d_solutions[3], c_solutions[3])
    @test evaluate_depths_new(Q[4], B[4], Ks[4], iF[4]) ≈ (d_solutions[4], c_solutions[4])
end

