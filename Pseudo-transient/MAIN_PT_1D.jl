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
    ρ      = 1000.0,   # kg/m^3
    useless = 1.0
)

# Physical parameters
Q = 20.0                # discharge [m^3/s]
B = 10.0                # channel width [m]
Ks = 30.0               # Strickler coefficient [m^(1/3)/s]
iF = 0.001              # channel slope [-]
L = 1000.0              # channel length [m]
n_points = 10           # number of points

# Numerical parameters
dx = L/(n_points-1)                 # space step [m]
xb = collect(dx/2:dx:L-dx/2)        # coordinate of the barycenters  (n_points)
xi = collect(0:dx:L)                # coordinate of the interfaces (n_points+1)
ρ_hat_bar = 1.0                     # parameter defining the pseudo-transient method

# Initialize arrays
b    = zeros(n_points)                # bed elevation [m] BARYCENTERS
h    = ones(n_points)                 # water depth [m]  BARYCENTERS (INITIAL VALUES ALREADY INITIALIZED)
q    = zeros(n_points+1)
u    = zeros(n_points+1)
iF_v = zeros(n_points+1)
dudx = zeros(n_points)

# Build the bed
build_bed!(b, x, iF, n_points)

# Evaluate the normal depth for the selected discharge
hL = evaluate_depth(Q, B, Ks, iF)
@infiltrate 

# set the initial conditions
h[1] = hL
q .= Q
# u .= Q./h

iter = 1; 
