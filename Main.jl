include("./QuantumIsing.jl")

using Plots
using ..QuantumIsing

hxs = 0:0.1:1.5
hsnew = 1.5:0.01:2
hs = 10 .^ range(-2., stop=1/2, length=50);
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

M1 = zeros(length(hs));
M2 = zeros(length(hs));
@time for (i,h) in enumerate(hs)
    vals1, vecs1 = generate_eigs(N=16, hx=0.5, hz=h)
    #vals2, vecs2 = generate_eigs(N=16, hx=h, hz=0.3)
    groundstate1 = @view vecs1[:,1]
    #groundstate2 = @view vecs2[:,1]
    M1[i] = magnetization(groundstate1,z=false)
    #M2[i] = magnetization(groundstate2)
    println(M1)
end
plot!(p, hs, M1, marker=:circle);
#plot!(p, hs, M2, marker=:square);

p

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

