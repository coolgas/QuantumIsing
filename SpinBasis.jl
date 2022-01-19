module SpinBasis
export Sx, Sy, Sz, Sx′, Sy′, Sz′, iS, iN, vg, ve
Sx = (1.0/2.0) * [0.0 1.0; 1.0 0.0]
Sy = (1.0/2.0) * [0.0 -1.0*1im; 1.0*1im 0.0]
Sz = (1.0/2.0) * [1.0 0.0; 0.0 -1.0]
iS = [1.0 0.0; 0.0 1.0]
iN = [1.0 0.0; 0.0 0.0]
vg = [0.0; 1.0]
ve = [1.0; 0.0]

# Rotated w.r.t y-axis π/2 degrees
Sx′ = (1.0/2.0) * [1.0 0.0;0.0 -1.0]
Sy′ = (1.0/2.0) * [0.0 -1.0*1im;1.0*1im 0.0]
Sz′ = (1.0/2.0) * [0.0 -1.0;-1.0 0.0]
end
