using Test, Printf, IJulia
include("../functions.jl")
include("../fun_opt.jl")

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

#@test evaluate_depths_new(Q[1], B[1], Ks[1], iF[1]) ≈ (d_solutions[1], c_solutions[1])
# d = c = 0.0
@testset "evaluate_depths_new" begin
    @test begin
    d, c = evaluate_depths_new(Q[1], B[1], Ks[1], iF[1]) 
    isapprox(d, d_solutions[1]; atol=1.e-5) && isapprox(c, c_solutions[1]; atol=1.e-5); end
    @test begin
    d, c = evaluate_depths_new(Q[2], B[2], Ks[2], iF[2])
    isapprox(d, d_solutions[2]; atol=1.e-5) && isapprox(c, c_solutions[2]; atol=1.e-5); end
    @test begin
    d, c = evaluate_depths_new(Q[3], B[3], Ks[3], iF[3])
    isapprox(d, d_solutions[3]; atol=1.e-5) && isapprox(c, c_solutions[3]; atol=1.e-5); end
    @test begin
    d, c = evaluate_depths_new(Q[4], B[4], Ks[4], iF[4])
    isapprox(d, d_solutions[4]; atol=1.e-5) && isapprox(c, c_solutions[4]; atol=1.e-5); end
end