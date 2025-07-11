using IJulia, Plots, Printf, Infiltrator 
include("functions.jl")

# ---------------------- PROBLEM DEFINITION -----------------------

n_reaches = 3               # number of reaches in the channel
L  = [1000.0, 2000.0, 1000.0]       # length of the channel [m]
Q  = 100.0                  # discharge [m^3/s] Q = Ω Ks R^2/3 iF^1/2
B  = [50.0, 50.0, 50.0]           # width of the channel [m]
Ks = [40.0, 40.0, 40.0]           # Strickler coefficient [m^(1/3)/s] Ks = 1/n
iF = [0.01, 0.001, 0.0005]          # slope of the channel
dx = 10.0                   # distance between the points of the channel [m]
g  = 9.81                   # gravity acceleration [m/s^2]

# -------------------------- GEOMETRY -----------------------------

n_points = Array{Int64}(undef, n_reaches)                                     # number of points in the channel
for n_rea in 1:n_reaches
    # Check the lengths of the channel are multiple of dx
    L[n_rea] % dx ≠ 0 && error("The length of the channel must be a multiple of dx")
    n_points[n_rea] = Int(L[n_rea]/dx) + 1
    @printf("Number of points along reach %d: %d \n", n_rea, n_points[n_rea])
end

n_points_total = sum(n_points)                              # total number of points in the channel
@printf("Total number of points in the channel: %d \n", n_points_total)

B_v = Array{Float64}(undef, n_points_total)                                  # temporary array for the widths
KS  = Array{Float64}(undef, n_points_total)                               # temporary array for the Strickler coefficients
IF  = Array{Float64}(undef, n_points_total)  

dX = dx .* ones(n_points_total-1)                            # vector of distances between the points
B_v[1:n_points[1]] = B[1] .* ones(n_points[1])                     # vector of widths, same length as X
KS[1:n_points[1]] = Ks[1] .* ones(n_points[1])                     # vector of Strickler coefficients
IF[1:n_points[1]] = iF[1] .* ones(n_points[1])                     # vector of slopes
for n_rea in 1:n_reaches-1
    location = sum(n_points[1:n_rea]) + 1                     # location of the first point of the next reach
    dX[location] = 0.0 # set the last distance to zero, to avoid problems in the integration
    B_v[sum(n_points[1:n_rea])+1: sum(n_points[1:n_rea+1])] = B[n_rea+1] .* ones(n_points[n_rea+1]) # vector of widths, same length as X
    KS[sum(n_points[1:n_rea])+1 : sum(n_points[1:n_rea+1])] = Ks[n_rea+1] .* ones(n_points[n_rea+1]) # vector of Strickler coefficients
    IF[sum(n_points[1:n_rea])+1 : sum(n_points[1:n_rea+1])] = iF[n_rea+1] .* ones(n_points[n_rea+1]) # vector of slopes
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



X = Array{Float64}(undef, n_points_total)                       # initialize the array for the x coordinates of the points
Z = Array{Float64}(undef, n_points_total)                       # initialize the arrays for the x and z coordinates of the points
for n_rea in 1:n_reaches
    x = LinRange(0, L[n_rea], n_points[n_rea])                  # x coordinates of the points
    z = zeros(n_points[n_rea])                                  # initialize z coordinates of the points
    build_bed!(z, x, dx, iF[n_rea], n_points[n_rea])            # build the bed of the channel
    if n_rea == 1
        X[1:n_points[n_rea]] = x .+ scale[n_rea, 1]                                                        # first reach, assign the x coordinates
        Z[1:n_points[n_rea]] = z .+ scale[n_rea, 2]                                        # first reach, assign the z coordinates with the scale
    else
        X[sum(n_points[1:n_rea-1])+1 : sum(n_points[1:n_rea])] = x .+ scale[n_rea, 1]   # append the x coordinates
        Z[sum(n_points[1:n_rea-1])+1 : sum(n_points[1:n_rea])] = z .+ scale[n_rea, 2]   # append the z coordinates with the scale
    end
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
    d_uniform[n_rea], d_critical[n_rea] = evaluate_depths(Q, B[n_rea], L[n_rea], Ks[n_rea], iF[n_rea])    # evaluate the uniform and critical depths
    # append the critical depth vector
    println("Reach ", n_rea, ": Uniform depth = ", d_uniform[n_rea], " m, Critical depth = ", d_critical[n_rea], " m")
end

# fill the critical depth vector with the values for each reach
d_critical_vector = zeros(n_points_total)                     # initialize the vector for the critical depths
d_critical_vector[1:n_points[1]] = d_critical[1] .* ones(n_points[1]) # first reach, assign the critical depth
for n_rea in 2:n_reaches
    d_critical_vector[sum(n_points[1:n_rea-1])+1 : sum(n_points[1:n_rea])] = d_critical[n_rea] .* ones(n_points[n_rea]) # append the critical depth
end
# d_critical_vector = n_rea == 1 ? d_critical[n_rea] .* ones(n_points[n_rea]) : vcat(d_critical_vector, d_critical[n_rea] .* ones(n_points[n_rea])) 


# The BC are always the uniform flow depth:
eL = d_uniform[begin] + Q^2 / (2 * g * (B[begin] * d_uniform[begin])^2) # water elevation at the left boundary [begin]
eR = d_uniform[end] + Q^2 / (2 * g * (B[end] * d_uniform[end])^2) # water elevation at the left boundary [end]
println("Left E BC = ", eL, " m, Right E BC = ", eR, " m")

# ------------------------- INTEGRATION ----------------------------- 
# NB: here we have to solve for E and immediately evaluate the depths from the energy


# Solve the energy equation from left to right (downward)
e_dw = zeros(n_points_total); e_dw[1] = eL                      # energy E at the points of the channel + IC
d_fast = zeros(n_points_total); d_fast[1] = d_uniform[begin]    # depth from the downward energy + IC
for n in 2:n_points_total
    if d_fast[n-1] == 0.0
        @printf("Warning: d_fast[%d] = %.2f m, set to critical: d_fast[%d] = %.3f m \nBreaking the supercritical computation  \n", n-1, d_fast[n-1], n-1, d_critical_vector[n-1])
        d_fast[n-1] = d_critical_vector[n-1] # set the depth to the critical depth
        break
    end
    e_dw[n] = e_dw[n-1] + dX[n-1] * dEdx(IF[n-1], Q, B_v[n-1], d_fast[n-1], KS[n-1])
    # Evaluate the depth from the energy
    d_fast[n] = E2d(d_critical_vector[n], Q, B_v[n], e_dw[n], style="supercritical")
    @printf("n = %d, e_dw[n] = %f, d_fast[n] = %.4f, X[n] = %.2f \n", n, e_dw[n], d_fast[n], X[n])
end
println("Supercritical computation: Left downward Energy = ", e_dw[1], " m, Right downward Energy = ", e_dw[n_points_total], " m")


# Solve the energy equation from right to left (upward)
e_uw = zeros(n_points_total); e_uw[end] = eR                    # energy E at the points of the channel + IC
d_slow = zeros(n_points_total); d_slow[end] = d_uniform[end]    # depth from the upward energy + IC
for n in n_points_total-1:-1:1
    if d_slow[n+1] == 0.0
        @printf("Warning: d_slow[%d] = %.2f m, set it to critical: d_slow[%d] = %.3f m \nbreaking the subcritical computation  \n", n+1, d_slow[n+1], n+1, d_critical_vector[n+1])
        d_slow[n+1] = d_critical_vector[n+1] # set the depth to the critical depth
        break
    end
    # Evaluate the Energy (forward explicit Euler)
    e_uw[n] = e_uw[n+1] - dX[n] * dEdx(IF[n+1], Q, B_v[n+1], d_slow[n+1], KS[n+1])
    @infiltrate false
    d_slow[n] = E2d(d_critical_vector[n], Q, B_v[n], e_uw[n], style="subcritical")
    @printf("n = %d, e_uw[n] = %f, d_slow[n] = %.4f, X[n] = %.2f \n", n, e_uw[n], d_slow[n], X[n])
end
println("Subcritical computation: Left upward Energy = ", e_uw[1], " m, Right upward Energy = ", e_uw[n_points_total], " m")

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

y_limits = (min(minimum(Z), minimum(e_dw)) - 0.05, max(maximum(Z), maximum(Z+e_dw)) + 0.05)
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
    grid=true, xlims=(900, 1050), ylims=(Z[100]-2, Z[100]+d_slow[100]+1), color="black", linewidth=2,
    size=(1400, 700))
p2 = plot!(X, Z + d_fast, label="Supercritical depths", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="red", linewidth=2)
p2 = plot!(X, Z + d_slow, label="Subcritical depths", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="blue", linewidth=2)
p2 = plot!(X, Z + d, label="Water depth", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="green", linewidth=2)
p2 = plot!(X, Z + d_critical_vector, label="Critical depth", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="gray", linewidth=2, linestyle=:dash)

display(plot(p1, p2, layout=(2,1), size=(1400, 800)))

