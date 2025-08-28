# This file is to test whether the Pseudo-Transient method works for simple hyperbolic eq.
using IJulia, Plots, Printf, LinearAlgebra, BenchmarkTools, Infiltrator, JLD2, SpecialFunctions

# functions
f(q) = @. 0.5 * q^2
function eig(q)
    ε = 1.e-7
    out = similar(q)
    out = @. (f(q+ε) - f(q-ε)) / (2ε)
    return out
end

@views function Lax_Friedrics!(q, dtdx, ql, qr)
    Nx   = length(q)
    qnew = zeros(Nx)
    fm   = zeros(Nx)
    fp   = zeros(Nx)
    for i ∈ 1:Nx
        if i == 1
            fm = 0.5*(f(q[i])   + f(ql)  ) - 0.5*dtdx*(q[i]   - ql  ) 
            fp = 0.5*(f(q[i+1]) + f(q[i])) - 0.5*dtdx*(q[i+1] - q[i])
        elseif i == Nx
            fm = 0.5*(f(q[i]) + f(q[i-1])) - 0.5*dtdx*(q[i] - q[i-1])
            fp = 0.5*(f(qr)   + f(q[i]  )) - 0.5*dtdx*(qr   - q[i])
        else # Lax-Friedrics flux
            fm = 0.5*(f(q[i]) + f(q[i-1])) - 0.5*dtdx*(q[i] - q[i-1])
            fp = 0.5*(f(q[i+1]) + f(q[i])) - 0.5*dtdx*(q[i+1] - q[i])
        end
        qnew[i] = q[i] - dtdx * (fp - fm)
    end
    q .= qnew
end

@views function Pseudotransient_1D!(q, dt, dx, ql)
    # Nomenclature:
    # q old = q n-1
    # q = q star
    # q new = q n

    Nx    = length(q)
    res   = zeros(Nx-1)
    qnew  = zeros(Nx)
    qold  = copy(q)
    ρ_hat = 1000.0
    iter  = 0; max_iter = 1000; tol = 1.e-8; err = 2*tol
    while err ≥ tol && iter ≤ max_iter
        # Automatically upwind if q > 0 (that's the case here)
        qnew[1] = q[1] - 1/ρ_hat/dt * (q[1] - qold[1]) - 1/ρ_hat/dx * (f(q[1]) - f(ql))
        for i ∈ 2:Nx
            qnew[i] = q[i] - 1/ρ_hat/dt * (q[i] - qold[i]) - 1/ρ_hat/dx * (f(q[i]) - f(q[i-1]))
        end

        # Convergence check
        if iter%25 == 0
            res .= abs.( 1/ρ_hat/dt * (qnew[2:end] - qold[2:end]) + 1/ρ_hat/dx * (f(q[2:end]) - f(q[1:end-1])) )
            err = maximum(res)
            @printf("iter = %d, err = %.5e \n", iter, err)     
        end

        iter += 1         
        q .= qnew
    end
end

@views function Burger_1D(solutor, saving, loading, IC; n_vis=5) 

    # Numerics
    Nx = 500
    xl, xr = -1.0, 1.0
    dx = (xr - xl) / Nx
    x = collect(xl+0.5*dx:dx:xr-0.5*dx)

    time = 0.0
    t_end = 0.5
    CFL = 0.9
    NMAX = 1000

    # Initial data
    q  = zeros(Nx)
    # fm = zeros(Nx)
    # fp = zeros(Nx)
    # q .= exp.(-0.5*x.^2 / 0.1^2) 
    if IC == "shock"
        q = [exp(-0.5*ix^2/0.1^2) for ix in x]
    elseif IC == "rarefaction"
        x0   = 0.0     # centro della transizione
        width = 0.1    # controlla la pendenza (più piccolo = più ripida)
        ymin = 0.1
        ymax = 1.0
        q = [ymin + (ymax - ymin) * 0.5 * (1 + erf((ix - x0)/width)) for ix in x]
    end
    ql = q[1]; qr = q[end]  # BCs

    # plot IC
    pltIC = plot(x, q, label="t = $time", xlabel="x", ylabel="q", title="1D Burger IC")
    # display(pltIC)

    @infiltrate false

    # Time loop
    for n ∈ 1:NMAX
        amax = maximum(abs.(eig(q)))
        dt = dx*CFL/amax
        if time+dt > t_end
            dt = t_end - time
        end
        time ≥ t_end && break

        if (solutor==1)
            
        elseif (solutor==2)
            dtdx = dt/dx
            Lax_Friedrics!(q, dtdx, ql, qr)
        elseif (solutor==2)
            Pseudotransient_1D!(q, dt, dx, ql)
        end

        # Update Time and solution
        time += dt

        if (n%n_vis == 0)
        plt = plot(x, q, xlabel="x", ylabel="q", title="1D Burger, t = $time")
        display(plt)
        end
    end

    # Load and plot results with reference solution if available
    if loading # with reference solution
        if IC == "shock"
            @load "t05_xq_shock.jld2" xx qq
        elseif IC == "rarefaction"
            @load "t05_xq_rarefaction.jld2" xx qq
        end
        plt = plot(xx, qq, xlabel="x", ylabel="q", title="1D Burger, t = $time", label="Lax-Friedrics")
        plt = plot!(x, q, xlabel="x", ylabel="q", label="Pseudo-transient")
        display(plt)
        savefig(plt, "Burger1D_rarefaction.png")
    else # without reference solution
        plt = plot(x, q, xlabel="x", ylabel="q", title="1D Burger, t = $time")
        display(plt)
    end

    # Save results
    if saving
        xx, qq = copy(x), copy(q)
        @save "t05_xq_shock.jld2" xx qq
    end
end

# Select the solutor for solving the equation
# 1 = Godononv
# 2 = Lax-Friedrichs
# 3 = Pseudotransient
solutor = 1
saving = false
loading = true
IC = "rarefaction" # "shock" or "rarefaction"
Burger_1D(solutor, saving, loading, IC, n_vis=1000)



