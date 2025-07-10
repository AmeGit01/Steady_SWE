# THERE COULD BE PROBLEMS WITH println

# ---------------------- PROBLEM DEFINITION -----------------------

L  = 1000.0       # length of the channel [m]
Q  = 100.0        # discharge [m^3/s] Q = Ω Ks R^2/3 iF^1/2
B  = 50.0         # width of the channel [m]
Ks = 40.0         # Strickler coefficient [m^(1/3)/s] Ks = 1/n
iF = 0.001         # slope of the channel
dx = 1.0          # distance between the points of the channel [m]
g  = 9.81         # gravity acceleration [m/s^2]

# -------------------------- GEOMETRY -----------------------------

# Check the length of the channel is a multiple of dx
L % dx ≠ 0 && error("The length of the channel must be a multiple of dx")
n_points = Int(L/dx) + 1                                    # number of points in the channel
println("Number of points along the channel: ",  n_points)

X = LinRange(0, L, n_points)                                # x coordinates of the points
Z = zeros(n_points)                                         # initialize z coordinates of the points
build_bed!(Z, X, dx, iF, n_points) # build the bed of the channel
println("first and last x coordinates: ", X[begin], " m, ", X[end], " m")
println("first and last z coordinates: ", Z[begin], " m, ", Z[end], " m")

# --------------------- BOUNDARY CONDITIONS ------------------------

# Evaluate the uniform and critical depths
d_uniform, d_critical = evaluate_depths(Q, B, L, Ks, iF)    # evaluate the uniform and critical depths
d_critical_vector = d_critical[1] * ones(length(X))         # vector of critical depths, same length as X
println("Uniform depths: ", d_uniform, " m, Critical depths: ", d_critical, " m")

# The BC are always the uniform flow depth:
eL = d_uniform + Q^2 / (2 * g * (B * d_uniform)^2) # water elevation at the left boundary [begin]
eR = d_uniform + Q^2 / (2 * g * (B * d_uniform)^2) # water elevation at the left boundary [end]

# ------------------------- INTEGRATION -----------------------------

# Solve the energy equation from left to right (downward)
e_dw = zeros(n_points)                          # energy E at the points of the channel
e_dw[1] = eL                                    # initial condition at the left boundary
for n in 2:n_points
    # Evaluate the hydraulic radius:
    Rh = (B[1] * d_uniform[1])/(2*d_uniform[1] + B[1])
    # Evaluate the Energy (forward explicit Euler)
    e_dw[n] = e_dw[n-1] - dx * (iF[1] - (Q/(Ks * B[1] * d_uniform[1]))^2 * Rh^(-4/3))
end
println("Left downward Energy = ", e_dw[1], " m, Right downward Energy = ", e_dw[n_points], " m")

# Solve the energy equation from right to left (upward)
e_uw = zeros(n_points)                          # energy E at the points of the channel
e_uw[n_points] = eR                              # initial condition at the right boundary        
for n in n_points-1:-1:1
    # Evaluate the hydraulic radius:
    Rh = (B[1] * d_uniform[1])/(2*d_uniform[1] + B[1])
    # Evaluate the Energy (forward explicit Euler)
    e_uw[n] = e_uw[n+1] - dx * (iF[1] - (Q/(Ks * B[1] * d_uniform[1]))^2 * Rh^(-4/3))
end
println("Left upward Energy = ", e_dw[1], " m, Right upward Energy = ", e_dw[n_points], " m")

# ------------------------- DEPTHS EVALUATION -----------------------

# From the solution of the energy abstract two depths:
d_fast = zeros(n_points)            # depth from the downward energy
d_slow = zeros(n_points)            # depth from the upward energy   
for n = 1:n_points
    d_fast[n] = E2d_single_reach(d_critical, Q, B, e_dw[n], style="supercritical")
    d_slow[n] = E2d_single_reach(d_critical, Q, B, e_uw[n], style="subcritical")
end

# ------------------------- SPINTA EVALUATION -----------------------

# evaluate the Spinta for both the supercritical and subcritical depths
Spinta_fast = zeros(n_points)                   # Spinta for the supercritical depths  
Spinta_slow = zeros(n_points)                   # Spinta for the subcritical depths
for n = 1:n_points
    Spinta_fast[n] = Spinta(Q, B, d_fast[n])    # Spinta for the supercritical depths
    Spinta_slow[n] = Spinta(Q, B, d_slow[n])    # Spinta for the subcritical depths
end
println("Spinta for the supercritical depths: ", Spinta_fast[1], " N/m, with depth: ", d_fast[1], " m")
println("Spinta for the subcritical depths: ", Spinta_slow[1], " N/m, with depth: ", d_slow[1], " m")

# Chose the depth with the maximum Spinta between the two:
d = zeros(n_points)                              # final depth to be plotted
for n = 1:n_points
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