include("./QuantumIsing.jl")

using Plots
using ..QuantumIsing

figpath = "./figs";

hxs = 0:0.1:1.5
hsnew = 1.5:0.01:2
hs = 10 .^ range(-2., stop=0.5, length=100);
hzs = 0.001:0.1:1.2;
Ns = 2:2:15;

p = plot();

@time for N in Ns
    M = zeros(length(hs))
    for (i,h) in enumerate(hs)
        vals, vecs = generate_eigs(N=N, hx=0.5, hz=h)
        groundstate = @view vecs[:,1]
        #M[i] = magnetization(groundstate)
        M[i] = staggered_magnetization(groundstate)
    end
    plot!(p, hs, M, marker=:circle, label="N = $N",
        xlab="h", ylab="M(h)")
    println(M)
end

p

vals, vecs = generate_eigs(N=14, hx=0.5, hz=1.5);
groundstate = @view vecs[:,1];
correlation(groundstate,j=2,step=2)


p = plot(xlab="hz",ylab="d^2Sx/dhz^2", xlim=(0.5, 2.5), ylim=(-2, 2));
M1 = zeros(length(hs));
M2 = zeros(length(hs));
M3 = zeros(length(hs));
M4 = zeros(length(hs));
@time for (i,h) in enumerate(hs)
    vals1, vecs1 = generate_eigs(N=8, hx=0.5, hz=h, rotated=false)
    vals2, vecs2 = generate_eigs(N=10, hx=0.5, hz=h, rotated=false)
    vals3, vecs3 = generate_eigs(N=12, hx=0.5, hz=h, rotated=false)
    vals4, vecs4 = generate_eigs(N=14, hx=0.5, hz=h, rotated=false)
    groundstate1 = @view vecs1[:,1]
    groundstate2 = @view vecs2[:,1] 
    groundstate3 = @view vecs3[:,1]
    groundstate4 = @view vecs4[:,1] 
    M1[i] = magnetization(groundstate1, z=false)
    M2[i] = magnetization(groundstate2, z=false)
    M3[i] = magnetization(groundstate3, z=false)
    M4[i] = magnetization(groundstate4, z=false)
    print(M4)
end

h1, deri1 = derivative(x=hs, y=M1);
h12, deri12 = derivative(x=h1, y=deri1);
h2, deri2 = derivative(x=hs, y=M2);
h22, deri22 = derivative(x=h2, y=deri2);
h3, deri3 = derivative(x=hs, y=M3);
h32, deri32 = derivative(x=h3, y=deri3);
h4, deri4 = derivative(x=hs, y=M4);
h42, deri42 = derivative(x=h4, y=deri4);
#plot!(p, h, deri, marker=:circle, label="1st order");
plot!(p, h12, deri12, marker=:circle, label="N=8");
plot!(p, h22, deri22, marker=:circle, label="N=10");
plot!(p, h32, deri32, marker=:circle, label="N=12");
plot!(p, h42, deri42, marker=:circle, label="N=14");
#plot!(p, hs, M1, marker=:circle);
#plot!(p, hs, M2, marker=:square);

p

savefig(p, "~/Desktop/Hamiltonian/figs/fig1.png")

p = plot(xlab="hx", ylab="hz", xlim=(0,1.2), legend=false);
annotate!((0.25, 0.5, "Antiferromagnetic"));
annotate!(1.00, 1.5, "Paramagnetic");
hxs = 0.06:0.01:1.2;
hzs = 0:0.2:2;
#hxs = 10 .^ range(-2., stop=0, length=10)
@time for hz in hzs
    for i in 1:1:length(hxs)-1
        m1, m2 = critical_line(N=10, hx=hxs[i], hz=hz) 
        m1prime, m2prime = critical_line(N=10, hx=hxs[i+1],hz=hz)
        if m2prime > m2
            hx = (hxs[i+1]+hxs[i])/2
            println(hx," ",hz)
            plot!(p, [hx], [hz], color=:blue, marker=:circle)
            break
        end
    end
end

p

