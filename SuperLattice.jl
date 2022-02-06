include("./SpinBasis.jl")
include("./Moperator.jl")

module SuperLattice

using LinearAlgebra
using ArnoldiMethod
using SparseArrays
using Arpack
using ..SpinBasis
using ..Moperator

export spectrum, global_ground_energy, generate_eigs

function spectrum(;v1::Real=0, v2::Real=0, v3::Real=0, Ω::Real=0, k::Real=0)
    ϵ_k=(1/8)*sqrt(((v1+v2)*sin(k)+v3*sin(2*k))^2+(8*Ω+(v1+v2)*cos(k)+v3*cos(2k))^2)
    return ϵ_k
end

"""
This function will give us the global ground energy at ABC sector
"""
function global_ground_energy(;v1::Real=0, v2::Real=0, v3::Real=0, Ω::Real=0, N::Int=0)
    @assert N % 2 == 0
    energy = 0
    for n in 1:N/2
        k = (2*n-1)*π/N
        energy += -1 * spectrum(v1=v1, v2=v2, v3=v3, Ω=Ω, k=k)
    end
    return energy
end

"""
Calculating the eigenvalue and eigenvector of the ground state of matrix x.
"""
function eigen_sparse(x)
    decomp, history = partialschur(x, nev=1, which=SR()) # only solve for the ground state
    λ, ϕ = partialeigen(decomp)
    return λ, ϕ
end

function eigen_sparse_multi(x)
    decomp, history = partialschur(x, nev=6, which=SR())
    λ, ϕ = partialeigen(decomp)
    return λ, ϕ
end

"""
This function will give us the eigenvalues and eigenvectors of the Quantum Ising
model with N sites and external magnetic fields hx and hz. 
"""
function generate_eigs(; N::Int, hx::Real=1, hz::Real=1, l::Real=1, a::Real=1, C6::Real=1, gstate=true)
    Ho = spzeros(2^N, 2^N)

    # First consider the interactions between 1-site
    for i in 1:N-1
        if i % 2 == 0
            vz = C6/(l-a)^6
        else
            vz = C6/a^6
        end
        Ho = Ho + vz*moperator(Sx, i, N)*moperator(Sx, i+1, N)
    end

    # Impose the PBC on the 1-site interactions
    Ho = Ho + (C6/(l-a)^6)*moperator(Sx, N, N)*moperator(Sx, 1, N)

    # Then we consider the 2-site interaction
    for i in 1:N-2
        vz = C6/l^6
        Ho = Ho + vz*moperator(Sx, i, N)*moperator(Sx, i+2, N)
    end

    # Impose the PBC on the 2-site interactions
    Ho = Ho + (C6/l^6)*moperator(Sx, N-1, N)*moperator(Sx, 1, N)
    Ho = Ho + (C6/l^6)*moperator(Sx, N, N)*moperator(Sx, 2, N)

    # Below consider the transverse field
    Hx = spzeros(2^N, 2^N)
    for i ∈ 1:N
        Hx = Hx + moperator(Sx, i, N)
    end
    
    # Below consider the longitudinal field
    Hz = spzeros(2^N, 2^N)
    for i ∈ 1:N
        Hz = Hz + moperator(Sz, i, N)
    end

    # Combining the pair interaction and transverse interaction
    Hxz = (1/4)*Ho + hx*Hx - hz*Hz
    
    λ = ϕ = 0
    if gstate
        λ, ϕ = eigen_sparse(Hxz) # this will only give us the ground state
    else
        λ, ϕ = eigen_sparse_multi(Hxz) # this will give us ground state and exicited states
    end

    return λ, ϕ
end

end


