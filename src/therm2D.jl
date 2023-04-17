const _SIGMA = 5.670374419E-8



function create_identity_matrix(n::T) where T<:Integer
    id_mat = zeros(T,n,n)
    for i = 1:n
        id_mat[i,i] = 1
    end
    return id_mat
end



function tempsolver(m, vfmat, temp, epsilon)
    # old algorithm used in matlab (slightly optimized)
    # prepare matrix for LGS
    # J --> outgoing radiation of element
    # G --> incoming radiation of element
    mat = ((epsilon .- 1) .* vfmat)
    # id_mat = Matrix{Int64}(I, m.no_elements, m.no_elements) # LinearAlgebra Pkg
    id_mat = create_identity_matrix(m.no_elements)
    mat .= mat .+ id_mat
    # solve LGS
    B = (temp[:,1].^4) .* _SIGMA .* epsilon
    J = mat \ B
    # Bestimmung der einfallenden Strahlung G
    G = zeros(Float64,m.no_elements,1)
    for i = 1:m.no_elements
        G[i,1] = sum(J[:,1] .* vfmat[i, :])
    end
    # Qp Radiation per element
    area = [m.elements[i].length for i = m.elem2par[1].first:m.elem2par[end].last]
    Qp = (J[:,1] - G[:,1]) .* area[:,1]
    Ga = G[:,1] .* area[:,1]
    return Qp, Ga
end



function set_bc_part!(m, bcv, parts, bcp)
    # give all face elements one temperature
    for i = parts
        e1 = m.elem2par[i].first
        e2 = m.elem2par[i].last
        bcv[e1:e2,1] .= bcp
    end
end
