include("./SuperLattice.jl")

using .SuperLattice
using Plots

g(x) = spectrum(v1=-10/1.9^6, v2=-10/2.1^6, v3=-10/4^6, Ω=0.04, k=x);
f(x) = -spectrum(v1=-10/1.9^6, v2=-10/2.1^6, v3=-10/4^6, Ω=0.04, k=x);

p = plot(legend=false, xlab="k", ylab="ϵ_k");
plot!(p, f, -π, π);
plot!(p, g, -π, π);
p

global_ground_energy(v1=-10/4^6, v2=-10/6^6, v3=-10/10^6, Ω=0.04, L=16)

generate_eigs(N=16, hx=0, hz=0.04, l=10, a=4, C6=-10)
