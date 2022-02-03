include("./SuperLattice.jl")

using .SuperLattice
using Plots

g(x) = spectrum(v1=-15/1^6, v2=-15/3^6, v3=-15/4^6, Ω=1, k=x);
f(x) = -spectrum(v1=-15/1^6, v2=-15/3^6, v3=-15/4^6, Ω=1, k=x);

p = plot(legend=false, xlab="k", ylab="ϵ_k");
plot!(p, f, -π, π);
plot!(p, g, -π, π);
p
