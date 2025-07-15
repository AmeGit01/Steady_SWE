using IJulia, Plots, Printf, Infiltrator 
include("fun_opt.jl")
GC.gc()
# ---------------------- PROBLEM DEFINITION -----------------------

n_reaches = 3               # number of reaches in the channel
L  = [1000.0, 2000.0, 1000.0]       # length of the channel [m]
Q  = 100.0                  # discharge [m^3/s] Q = Ω Ks R^2/3 iF^1/2
B  = [50.0, 50.0, 50.0]           # width of the channel [m]
Ks = [40.0, 40.0, 40.0]           # Strickler coefficient [m^(1/3)/s] Ks = 1/n
iF = [0.05, 0.0005, 0.01]          # slope of the channel
dx = 1.0                   # distance between the points of the channel [m]
g  = 9.81                   # gravity acceleration [m/s^2]

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

B_v = Float64[]                                        # vector of widths
# KS  = Float64[]                                        # vector of Strickler coefficients
# IF  = Float64[]                                        # vector of slopes
for n_rea in 1:n_reaches
    append!(B_v, B[n_rea].*ones(n_points[n_rea]))      # vector of widths
    # append!(KS, Ks[n_rea].*ones(n_points[n_rea]))      # vector of Strickler coefficients
    # append!(IF, iF[n_rea].*ones(n_points[n_rea]))      # vector of slopes
end

println("size B_v =  ", size(B_v1), ", size KS = ", size(KS1), ", size IF = ", size(IF1), ", size dX = ", size(dX))

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

X = Float64[]
Z = Float64[]
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
    color="black", linewidth=2, size=(1000, 400)); # display(pbed)


px   = plot(X, label="X coordinates", xlabel="Point index", ylabel="X [m]",
    title="X coordinates of the points", grid=true, color="blue",
    xlims=(1, n_points_total), size=(1000, 400))
# display(plot(pbed, px, layout=(2,1), size=(1000, 800), legend=:topright)) 


# --------------------- BOUNDARY CONDITIONS ------------------------

# Evaluate the uniform and critical depths
d_uniform, d_critical = zeros(n_reaches), zeros(n_reaches)    # initialize the vectors for the uniform and critical depths
for n_rea in 1:n_reaches
    d_uniform[n_rea], d_critical[n_rea] = evaluate_depths(Q, B[n_rea], Ks[n_rea], iF[n_rea])    # evaluate the uniform and critical depths
    # append the critical depth vector
    println("Reach ", n_rea, ": Uniform depth = ", d_uniform[n_rea], " m, Critical depth = ", d_critical[n_rea], " m")
end

# check wether the reaches are subcritical or supercritical
subcritical = zeros(n_reaches) # initialize the vector for the subcritical reaches
for n_rea in 1:n_reaches
    if d_uniform[n_rea] > d_critical[n_rea]
        subcritical[n_rea] = 1 # subcritical reach
        println("Reach ", n_rea, " is subcritical")
    else
        subcritical[n_rea] = 0 # supercritical reach
        println("Reach ", n_rea, " is supercritical")
    end
end

# fill the critical depth vector with the values for each reach
d_critical_vector = Float64[]
for n_rea in 1:n_reaches
    append!(d_critical_vector, d_critical[n_rea] .* ones(n_points[n_rea]))
end

# evaluate the BC for each reach (considered separately), so that to impose them internally where requested
EL = zeros(n_reaches); dL = zeros(n_reaches) # left boundary condition for the energy
ER = zeros(n_reaches); dR = zeros(n_reaches) # right boundary condition for the energy
EL[1]         = Energy(Q, B[1],         d_uniform[1],         g); 
dL[1]         = d_uniform[1]  # water elevation at the left boundary [m]
ER[n_reaches] = Energy(Q, B[n_reaches], d_uniform[n_reaches], g); 
dR[n_reaches] = d_uniform[n_reaches] # water elevation at the right boundary [m]
for n_rea in 2:n_reaches
    discriminant = subcritical[n_rea] - subcritical[n_rea-1] # check if the reach is subcritical or supercritical
    if discriminant == 0 || discriminant == +1
        EL[n_rea]   = Energy(Q, B[n_rea], d_uniform[n_rea], g);     dL[n_rea]   = d_uniform[n_rea] # if the reach is subcritical, use the uniform depth
        ER[n_rea-1] = Energy(Q, B[n_rea-1], d_uniform[n_rea-1], g); dR[n_rea-1] = d_uniform[n_rea]
    elseif discriminant == -1
        EL[n_rea]   = Energy(Q, B[n_rea], d_critical[n_rea], g);     dL[n_rea]   = d_critical[n_rea] # if the reach is supercritical, use the critical depth
        ER[n_rea-1] = Energy(Q, B[n_rea-1], d_critical[n_rea-1], g); dR[n_rea-1] = d_critical[n_rea-1] # if the previous reach is supercritical, use the critical depth
    else
        error("Discriminant must be either 0, -1, or 1")
    end
end

# Print the BC for checking purposes
for n_rea in 1:n_reaches
    @printf("Reach %d: Left E BC = %.4f m, Right E BC = %.4f m \n", n_rea, EL[n_rea], ER[n_rea])
end 


# ------------------------- INTEGRATION ----------------------------- 
# NB: here we have to solve for E and immediately evaluate the depths from the energy

# Solve the energy equation from left to right (downward)
# here we have to solve the supercritical reaches, starting from the left boundary and
# going on until reaching the critical state
e_dw   = zeros(n_points_total)
d_fast = zeros(n_points_total) # depth from the downward energy
for n_rea in 1:n_reaches
    if subcritical[n_rea] == 1
        @printf("Reach %d is subcritical, skipping the supercritical computation \n", n_rea)
        continue # skip the supercritical computation for subcritical reaches
    else
        @printf("Supercritical computation for reach %d \n", n_rea)
    end
    n = n_rea==1 ? 1 : sum(n_points[1:n_rea-1])+1 # starting point for the reach

    # skip the computation if the reach is supercritical and the previous one was supercritical too
    if n_rea>1 && (subcritical[n_rea-1] == 0 && subcritical[n_rea] == 0); continue; end 
    e_dw[n] = EL[n_rea]
    d_fast[n] = dL[n_rea]
    @printf("Imposed BC: depth = %.4f m,  energy = %.4f m \n", d_fast[n], e_dw[n])

    while n < n_points_total
        n += 1
        e_dw[n] = e_dw[n-1] + dX[n-1] * dEdx(IF[n-1], Q, B_v[n-1], d_fast[n-1], KS[n-1])
        d_fast[n] = analytical_E2d(d_critical_vector[n], e_dw[n], Q, B_v[n], g, n, style="supercritical")

        # @printf("n = %d, e_dw[n] = %f, d_fast[n] = %.4f, X[n] = %.2f \n", n, e_dw[n], d_fast[n], X[n])
    end
end
println("Supercritical computation: Left downward Energy = ", e_dw[1], " m, Right downward Energy = ", e_dw[n_points_total], " m")

# Solve the energy equation from right to left (upward)
# here we have to solve the subcritical reaches, starting from the right boundary and
# going on until reaching the critical state
e_uw   = zeros(n_points_total) 
d_slow = zeros(n_points_total) 
for n_rea in n_reaches:-1:1
    if subcritical[n_rea] == 0
        @printf("Reach %d is supercritical, skipping the subcritical computation \n", n_rea)
        continue # skip the supercritical computation for subcritical reaches
    else
        @printf("Subcritical computation for reach %d \n", n_rea)
    end
    n = n_rea==1 ? n_points[n_rea] : sum(n_points[1:n_rea])
    if n_rea < n_reaches && (subcritical[n_rea+1] == 1 && subcritical[n_rea] == 1); continue; end 
    e_uw[n] = ER[n_rea]
    d_slow[n] = dR[n_rea] # initial guess for the depth
    while n > 1
        n -= 1
        e_uw[n] = e_uw[n+1] - dX[n] * dEdx(IF[n+1], Q, B_v[n+1], d_slow[n+1], KS[n+1])
        d_slow[n] = analytical_E2d(d_critical_vector[n], e_uw[n], Q, B_v[n], g, n, style="subcritical")

        # @printf("n = %d, e_uw[n] = %f, d_slow[n] = %.4f, X[n] = %.2f \n", n, e_uw[n], d_slow[n], X[n])
    end
end

# ------------------------- SPINTA EVALUATION -----------------------

# evaluate the Spinta for both the supercritical and subcritical depths
Spinta_fast = zeros(n_points_total)                   # Spinta for the supercritical depths  
Spinta_slow = zeros(n_points_total)                   # Spinta for the subcritical depths
d = zeros(n_points_total)                             # final depth to be plotted
for n = 1:n_points_total
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
p1 = plot!(X, Z + d_critical_vector, label="Critical depth", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="gray", linewidth=2, linestyle=:dash)

p2 = plot(X, Z, label="Bed", xlabel="x [m]", ylabel="z [m]", title="Bed of the channel", legend=:topright,
    grid=true, xlims=(900, 1050), ylims=(Z[1000]-1, Z[1000]+d_slow[100]+4), color="black", linewidth=2,
    size=(1400, 700))
p2 = plot!(X, Z + d_fast, label="Supercritical depths", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="red", linewidth=2)
p2 = plot!(X, Z + d_slow, label="Subcritical depths", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="blue", linewidth=2)
p2 = plot!(X, Z + d, label="Water depth", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="green", linewidth=2)
p2 = plot!(X, Z + d_critical_vector, label="Critical depth", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="gray", linewidth=2, linestyle=:dash)

# display(plot(p1, p2, layout=(2,1), size=(1400, 800)))

