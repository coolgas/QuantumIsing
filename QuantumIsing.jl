include("./SpinBasis.jl")
include("./Moperator.jl")

module QuantumIsing 

using ArnoldiMethod
using SparseArrays
using Arpack
using ..SpinBasis
using ..Moperator

export generate_eigs, magnetization, critical_line, staggered_magnetization, derivative, correlation

"""
Generate the a base of the Hamiltonian in terms of the binary numbers 0 and 1.
"""
bit_rep(num::Integer, N::Integer) = BitArray(parse(Bool, i) for i in string(num, base=2, pad=N))

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

# Taking the derivative
function derivative(;x=Array{Float64}(undef), y=Array{Float64}(undef), step=1)
    @assert length(x) == length(y)
    @assert step <= length(y)-1
    X = [x[i] for i in 1:step:length(x)]
    Y = [y[i] for i in 1:step:length(y)]
    derivative = [(Y[i+1]-Y[i])/(X[i+1]-X[i]) for i in 1:1:length(X)-1]
    return X[1:length(X)-1], derivative
end
        
"""
This function will give us the eigenvalues and eigenvectors of the Quantum Ising
model with N sites and external magnetic fields hx and hz. 
"""
function generate_eigs(; N::Int, hx=1, hz=1, gstate=true, rotated=false)
    Ho = spzeros(2^N, 2^N)
    # consider the rotated case first
    if rotated
        for i ∈ 1:N
            for j ∈ 1:N
                if abs(i-j) == 1
                    vz = 1
                elseif abs(i-j) == N-1 # apply the PBC
                    vz = 1
                else
                    vz = 0
                end
                Ho = Ho + vz * moperator(Sz′,i,N) * moperator(Sz′,j,N)
            end
        end
        # Below consider the transverse field
        Hx = spzeros(2^N, 2^N)
        for i ∈ 1:N
            Hx = Hx + moperator(Sx′, i, N)
        end
        
        # Below consider the longitudinal field
        Hz = spzeros(2^N, 2^N)
        for i ∈ 1:N
            Hz = Hz + moperator(Sz′, i, N)
        end

    # Below considers the case of pair interaction within distance
    else
        for i ∈ 1:N
            for j ∈ 1:N
                if abs(i-j) == 1
                    vz = 1
                elseif abs(i-j) == N-1 # apply the PBC
                    vz = 1
                else
                    vz = 0
                end
                Ho = Ho + vz * moperator(Sz,i,N) * moperator(Sz,j,N)
            end
        end
        
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
    end
    # Combining the pair interaction and transverse interaction
    Hxz = Ho - hx*Hx - hz*Hz
    
    λ = ϕ = 0
    if gstate
        λ, ϕ = eigen_sparse(Hxz) # this will only give us the ground state
    else
        λ, ϕ = eigen_sparse_multi(Hxz) # this will give us ground state and exicited states
    end

    return λ, ϕ
end

# Calculate the correlation functions at site j with custmized steps
function correlation(state; j::Int, step::Int, z=true)
    N = Int(log2(length(state)))
    @assert 0 < j <= N
    @assert 0 <= step < N
    @assert 0 < j + step <= N
    corr = 0 # the output correlation
    if z
        corr = state' * moperator(Sz,j,N) * moperator(Sz,j+step,N) * state
    else
        corr = state' * moperator(Sx,j,N) * moperator(Sx,j+step,N) * state
    end
    
    return corr
end


# This will reproduce the phase diagram w.r.t hz and hx.
function critical_line(;N::Int, hx::Real, hz::Real, rotated=false)
    λ, ϕ = generate_eigs(N=N, hx=hx, hz=hz, gstate=false, rotated=rotated)
    c = collect(zip(sortperm(λ), λ[sortperm(λ)]))

    ϕ0 = ϕ[:,c[1][1]]
    ϕ1 = ϕ[:,c[2][1]]
    ϕ2 = ϕ[:,c[3][1]]
    
    Ho = spzeros(2^N, 2^N)

    if rotated
        for i ∈ 1:N
            for j ∈ 1:N
                if abs(i-j) == 1
                    vz = 1
                elseif abs(i-j) == N-1 # apply the PBC
                    vz = 1
                else
                    vz = 0
                end
                Ho = Ho + vz * moperator(Sz′,i,N) * moperator(Sz′,j,N)
            end
        end

        Hx = spzeros(2^N, 2^N)
        for i ∈ 1:N
            Hx = Hx + moperator(Sx′, i, N)
        end
        
        Hz = spzeros(2^N, 2^N)
        for i ∈ 1:N
            Hz = Hz + moperator(Sz′, i, N)
        end
 
    else
        for i ∈ 1:N
            for j ∈ 1:N
                if abs(i-j) == 1
                    vz = 1
                elseif abs(i-j) == N-1 # apply the PBC
                    vz = 1
                else
                    vz = 0
                end
                Ho = Ho + vz * moperator(Sz,i,N) * moperator(Sz,j,N)
            end
        end
        Hx = spzeros(2^N, 2^N)
        for i ∈ 1:N
            Hx = Hx + moperator(Sx, i, N)
        end
        
        Hz = spzeros(2^N, 2^N)
        for i ∈ 1:N
            Hz = Hz + moperator(Sz, i, N)
        end
    end

    m1 = hx*(ϕ0'*Hx*ϕ0-ϕ1'*Hx*ϕ1)+hz*(ϕ0'*Hz*ϕ0-ϕ1'*Hz*ϕ1)-ϕ0'*Ho*ϕ0+ϕ1'*Ho*ϕ1
    m2 = hx*(ϕ0'*Hx*ϕ0-ϕ2'*Hx*ϕ2)+hz*(ϕ0'*Hz*ϕ0-ϕ2'*Hz*ϕ2)-ϕ0'*Ho*ϕ0+ϕ2'*Ho*ϕ2

    return m1, m2
end

"""
This will calculate the magnetization of the state.
"""
function magnetization(state; z=true, i=1)
    N = Int(log2(length(state)))
    M = 0
    if z
        for i in 1:length(state)
            bstate = bit_rep(i-1, N)
            bstate_M = 0
            for spin in bstate
                bstate_M += (state[i]^2 * (spin ? 1 : -1))/N
            end
            #@assert abs(bstate_M) <= 1
            M += abs(bstate_M)
            #M += bstate_M

        end
    else
        M = state' * moperator(Sx,i,N) * state
    end

    return M
end

function staggered_magnetization(state)
    N = Int(log2(length(state)))
    SM = 0
    for i in 1:length(state)
        bstate = bit_rep(i-1,N)
        bstate_SM = 0
        for spin in bstate
            bstate_SM += (state[i]^2 * (-1)^(i) * (spin ? 1 : -1))/N
        end
        #@assert abs(bstate_SM) <= 1
        SM += abs(bstate_SM)
    end
    return SM
end

end
