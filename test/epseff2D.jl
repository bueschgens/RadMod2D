#= ######################################################
includes
###################################################### =# 

include("./models2D.jl")

#= ######################################################
inside functions
###################################################### =# 

function get_epsilon_effective(m, mat, p_hole; epsilon_ink = 0.1)
    # calculate epsilon effective
    # hole is part with epsilon effective / encl is enclosure
    p_encl = collect(1:m.npar)
    filter!(x->x != p_hole, p_encl)
    temp_encl = 0
    temp_hole = 1000 # Kelvin
    temp = zeros(m.nelem,1)
    set_bc_part!(m, temp, p_hole, temp_hole)
    set_bc_part!(m, temp, p_encl, temp_encl)
    epsilon_hole = 1.0
    n_epsilon_encl = round(Integer,1.0/epsilon_ink)
    epsilon_encl = Matrix{Float64}(undef,n_epsilon_encl,2)
    epsilon_encl[:,1] = [epsilon_ink*i for i =1:n_epsilon_encl]
    epsilon = zeros(m.nelem,1)
    set_bc_part!(m, epsilon, p_hole, epsilon_hole)
    for i = 1:n_epsilon_encl
        set_bc_part!(m, epsilon, p_encl, epsilon_encl[i,1])
        # calculation of Qp
        Qp, Ga = tempsolver(m, mat, temp, epsilon)
        area_hole = get_area_of_part(m, p_hole)
        G_hole = sum(Ga[m.elem2par[p_hole].first:m.elem2par[p_hole].last])
        Qp_hole = sum(Qp[m.elem2par[p_hole].first:m.elem2par[p_hole].last])
        # calculation of epsilon effective
        epsilon_eff = 1.0 - (G_hole / (area_hole * sigma * (temp_hole^4 - temp_encl^4)))
        epsilon_encl[i,2] = epsilon_eff
    end
    return epsilon_encl
end