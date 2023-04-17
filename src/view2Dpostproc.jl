"""
    compact_vfmat_to_parts(m::ConstModel, vfmat::Matrix{T}; normit::Bool = false) where T <: AbstractFloat


Compact the viewfactor matrix `mat` n_elements x n_elements to the parts of the model
`m` n_parts x n_parts. 

# Arguments
- `m::ConstModel`: model with parts
- `mat::Matrix{T}`: viewfactor matrix

# Keyword Arguments
- `normit::Bool`: normalize the viewfactor matrix to the sum of viewfactors of each element

# Returns
- `matp::Matrix{T}`: viewfactor matrix of parts
"""
function compact_vfmat_to_parts(m::ConstModel, vfmat::Matrix{T}; normit::Bool = false) where T <: AbstractFloat

    vfmatp = zeros(T, m.no_parts, m.no_parts)
    areap = [get_area_of_part(m, i) for i = 1:m.no_parts]

    for p1 = 1:m.no_parts, p2 = 1:m.no_parts
        vfsplit = vfmat[m.elem2par[p1].first:m.elem2par[p1].last, m.elem2par[p2].first:m.elem2par[p2].last]
        acol = [m.elements[i].length for i = m.elem2par[p1].first:m.elem2par[p1].last]
        asplit = repeat(acol, 1, (m.elem2par[p2].last-m.elem2par[p2].first+1))
        vfsplit = vfsplit .* asplit
        vfmatp[p1,p2] = sum(vfsplit) / areap[p1]
    end

    if normit
        controlp = sum(vfmatp, dims=2)
        vfmatp .= vfmatp ./ controlp
    end

    return vfmatp
end



"""
    save_vfmat_part_to_csv(mat, path, name; control = false, partmatrix=true))

write the viewfactor matrix `mat` n_parts x n_parts to a csv file in `path` with name
`name`.

# Arguments
- `mat::Matrix{Float64}`: viewfactor matrix

# Keyword Arguments
- `filepath::String`: path to save the file
- `filename::String`: name of the file
- `control::Bool`: add a control column to the csv file
- `partmatrix::Bool`: if true, the part names are added to the header of the csv file
  (matrix must be size n_parts x n_parts)
"""
function save_vfmat_part_to_csv(vfmat::Matrix{T}; filepath::String = ".",
        filename::String = "vfmat.csv", control = false, partmatrix=false,
        model = nothing)::Nothing where T <: AbstractFloat

    header = []

    if partmatrix
        header = Matrix{String}(undef, 1, size(vfmat,1))

        for i = 1:size(vfmat,1)
            header[1,i] = "p"*string(i)

            if !isnothing(model)
                header[1,i] *= ":"*model.part_names[i]
            end
        end
    end

    if control == true
        control = sum(vfmat, dims=2)
        vfmat = hcat(vfmat, control)
        header = hcat(header, "control")
    end

    fp = joinpath(filepath, filename)
    open(fp, "w") do io
        writedlm(io, header)
        writedlm(io, vfmat)
    end
end
