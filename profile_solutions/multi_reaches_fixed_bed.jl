using IJulia, Plots, Printf, Infiltrator 

# ---------------------- PROBLEM DEFINITION -----------------------

n_reaches = 2               # number of reaches in the channel
L  = [1000.0, 2000.0]       # length of the channel [m]
Q  = 100.0                  # discharge [m^3/s] Q = Ω Ks R^2/3 iF^1/2
B  = [50.0, 50.0]           # width of the channel [m]
Ks = [40.0, 40.0]           # Strickler coefficient [m^(1/3)/s] Ks = 1/n
iF = [0.01, 0.001]          # slope of the channel
dx = 1.0                    # distance between the points of the channel [m]
g  = 9.81                   # gravity acceleration [m/s^2]

# -------------------------- GEOMETRY -----------------------------

# Check the lengths of the channel are multiple of dx

n_points = Array{Int64}(undef, n_reaches)                                     # number of points in the channel
B_v, KS, IF, dX = Float64[], Float64[], Float64[], Float64                # vector of widths, same length as X
for n_rea in 1:n_reaches
    L[n_rea] % dx ≠ 0 && error("The length of the channel must be a multiple of dx")
    n_points[n_rea] = Int(L[n_rea]/dx) + 1
    println("Number of points along reach ", n_rea, ": ", n_points[n_rea])
    B_v = n_rea == 1 ? B[n_rea] .* ones(n_points[n_rea]) : vcat(B_v, B[n_rea] .* ones(n_points[n_rea])) # vector of widths, same length as X
    KS = n_rea == 1 ? Ks[n_rea] .* ones(n_points[n_rea]) : vcat(KS, Ks[n_rea] .* ones(n_points[n_rea])) # vector of Strickler coefficients
    IF = n_rea == 1 ? iF[n_rea] .* ones(n_points[n_rea]) : vcat(IF, iF[n_rea] .* ones(n_points[n_rea])) # vector of slopes
    dX = n_rea == 1 ? dx .* ones(n_points[n_rea]-1) : vcat(dX, dx .* ones(n_points[n_rea])) # vector of distances between the points
end
println("size B_v = ", size(B_v), ", size KS = ", size(KS), ", size IF = ", size(IF), ", size dX = ", size(dX))

for n_rea in 1:n_reaches-1
    location = sum(n_points[1:n_rea]) + 1                     # location of the first point of the next reach
    dX[location] = 0.0 # set the last distance to zero, to avoid problems in the integration
end
n_points_total = sum(n_points)                              # total number of points in the channel
println("Total number of points in the channel: ", n_points_total)

# Check the dX vector is corrected (must be present some zeros)
pdx = scatter(dX, label="dX", xlabel="Point index", ylabel="dX [m]",
    title="Distance between the points", grid=true, color="black",
    xlims=(1, n_points_total), size=(1000, 400))
display(pdx)    

# auxiliary vector
scale = zeros(n_reaches)                                    # scale vector for the z coordinates
for n_rea in n_reaches-1:-1:1
    scale[n_rea] = scale[n_rea+1] + L[n_rea] * iF[n_rea]                    # scale vector for the z coordinates
end
println("Scale vector: ", scale)

X = Float64[]; Z = Float64[]                                # initialize the arrays for the x and z coordinates of the points
for n_rea in 1:n_reaches
    x = LinRange(0, L[n_rea], n_points[n_rea])                     # x coordinates of the points
    z = zeros(n_points[n_rea])                              # initialize z coordinates of the points
    build_bed!(z, x, dx, iF[n_rea], n_points[n_rea])        # build the bed of the channel
    X = n_rea == 1 ? x : vcat(X, x)                             # append the x coordinates as it is
    Z = n_rea == 1 ? z .+ scale[n_rea] : vcat(Z, z .+ scale[n_rea])             # append the z coordinates with the scale
end

println("first and last x coordinates: ", X[1], " m, ", X[end], " m")
println("first and last z coordinates: ", Z[1], " m, ", Z[end], " m")

# --------------------- BOUNDARY CONDITIONS ------------------------

# Evaluate the uniform and critical depths
d_uniform, d_critical = zeros(n_reaches), zeros(n_reaches)    # initialize the vectors for the uniform and critical depths
for n_rea in 1:n_reaches
    d_uniform[n_rea], d_critical[n_rea] = evaluate_depths(Q, B[n_rea], L[n_rea], Ks[n_rea], iF[n_rea])    # evaluate the uniform and critical depths
    # append the critical depth vector
    d_critical_vector = n_rea == 1 ? d_critical[n_rea] .* ones(n_points[n_rea]) : vcat(d_critical_vector, d_critical[n_rea] .* ones(n_points[n_rea])) 
    println("Reach ", n_rea, ": Uniform depth = ", d_uniform[n_rea], " m, Critical depth = ", d_critical[n_rea], " m")
end

# The BC are always the uniform flow depth:
eL = d_uniform[begin] + Q^2 / (2 * g * (B[begin] * d_uniform[begin])^2) # water elevation at the left boundary [begin]
eR = d_uniform[end] + Q^2 / (2 * g * (B[end] * d_uniform[end])^2) # water elevation at the left boundary [end]
println("Left E BC = ", eL, " m, Right E BC = ", eR, " m")

# ------------------------- INTEGRATION ----------------------------- HERE THERE IS A PROBLEM
# NB: here we have to solve for E and immediately evaluate the depths from the energy

# Solve the energy equation from left to right (downward)
e_dw = zeros(n_points_total); e_dw[1] = eL                      # energy E at the points of the channel + IC
d_fast = zeros(n_points_total); d_fast[1] = d_uniform[begin]    # depth from the downward energy + IC
for n in 2:n_points_total
    # Evaluate the hydraulic radius:
    Rh = (B_v[n-1] * d_fast[n-1])/(2*d_fast[n-1] + B_v[n-1])
    # Evaluate the Energy (forward explicit Euler)
    temp = - dX[n-1] * (IF[n-1] - (Q/(KS[n-1] * B_v[n-1] * d_fast[n-1]))^2 * Rh^(-4/3))
    e_dw[n] = e_dw[n-1] + temp
    d_fast[n] = E2d(d_critical_vector[n], Q, B_v[n], e_dw[n], style="supercritical")
end
println("Supercritical computation: Left downward Energy = ", e_dw[1], " m, Right downward Energy = ", e_dw[n_points_total], " m")


# Solve the energy equation from right to left (upward)
e_uw = zeros(n_points_total); e_uw[end] = eR                    # energy E at the points of the channel + IC
d_slow = zeros(n_points_total); d_slow[end] = d_uniform[end]    # depth from the upward energy + IC
for n in n_points_total-1:-1:1
    # Evaluate the hydraulic radius:
    Rh = (B_v[n+1] * d_slow[n+1])/(2*d_slow[n+1] + B_v[n+1])
    # Evaluate the Energy (forward explicit Euler)
    temp = - dX[n] * (IF[n+1] - (Q/(KS[n+1] * B_v[n+1] * d_slow[n+1]))^2 * Rh^(-4/3))
    e_uw[n] = e_uw[n+1] + temp
    d_slow[n] = E2d(d_critical_vector[n], Q, B_v[n], e_uw[n], style="subcritical")
    println("n = ", n)
end
println("Subcritical computation: Left upward Energy = ", e_dw[1], " m, Right upward Energy = ", e_dw[n_points_total], " m")

#=
# ------------------------- SPINTA EVALUATION -----------------------

# evaluate the Spinta for both the supercritical and subcritical depths
Spinta_fast = zeros(n_points_total)                   # Spinta for the supercritical depths  
Spinta_slow = zeros(n_points_total)                   # Spinta for the subcritical depths
for n = 1:n_points_total
    Spinta_fast[n] = Spinta(Q, B_v[n], d_fast[n])    # Spinta for the supercritical depths
    Spinta_slow[n] = Spinta(Q, B_v[n], d_slow[n])    # Spinta for the subcritical depths
end
#println("Spinta for the supercritical depths: ", Spinta_fast[1], " N/m, with depth: ", d_fast[1], " m")
#println("Spinta for the subcritical depths: ", Spinta_slow[1], " N/m, with depth: ", d_slow[1], " m")

# Chose the depth with the maximum Spinta between the two:
d = zeros(n_points_total)                              # final depth to be plotted
for n = 1:n_points_total
    # use the supercritical depth if Spinta is greater, otherwise use the subcritical depth
    Spinta_fast[n] > Spinta_slow[n] ? d[n] = d_fast[n] : d[n] = d_slow[n] 
end

# --------------------------- PLOTS -----------------------------

y_limits = (min(minimum(Z), minimum(e_dw)) - 0.05, max(maximum(Z), maximum(Z+e_dw)) + 0.05)
p1 = plot(X, Z, label="Bed", xlabel="x [m]", ylabel="z [m]", title="Bed of the channel", legend=:topright,
    grid=true, xlims=(X[begin]-10, X[end]+10), ylims=y_limits, color="black", linewidth=2,
    size=(1400, 400))
p1 = plot!(X, Z + d_fast, label="Supercritical depths", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="red", linewidth=2)
p1 = plot!(X, Z + d_slow, label="Subcritical depths", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="blue", linewidth=2)
p1 = plot!(X, Z + d, label="Water depth", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="green", linewidth=2)
p1 = plot!(X, Z + d_critical_vector, label="Critical depth", xlabel="x [m]", ylabel="z [m]",
    title="Water elevation in the channel", grid=true, color="gray", linewidth=2, linestyle=:dash)
display(p1)
=#