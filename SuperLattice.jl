module SuperLattice

export spectrum

function spectrum(;v1::Real=0, v2::Real=0, v3::Real=0, Ω::Real=0, k::Real=0)
    ϵ_k=(1/8)*sqrt(((v1+v2)*sin(k)+v3*sin(2*k))^2+(8*Ω+(v1+v2)*cos(k)+v3*cos(2k))^2)
    return ϵ_k
end

end
