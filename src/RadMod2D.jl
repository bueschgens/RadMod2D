module RadMod2D

    # julia packages
    using DelimitedFiles
    using StaticArrays
    
    include("./struct2D.jl")
    export Point2D, Index2D

    include("./mesh2D.jl")
    export ModelRaw, Model
    export create_empty_model
    export add!, offset_model!
    export make_model_immutable
    export edge, rectangle, circle, circle_open, cosinus
    export element_analysis

    include("./raycast2D.jl")
    export create_randomly_occupied_tiles
    export get_max_steps
    export tilewalk_with_check_for_occ_tiles
    export tilewalk_with_return, tilewalk_with_return!

    include("./view2D.jl")
    export existing_vf!
    export blocking_vf!
    export create_tiles, check_tile_occupation
    export blocking_vf_with_tiles!
    export calculating_vf!
    export get_area_of_part
    export vfmat_to_parts
    # 2 elems
    export existing_vf_2elem
    export blocking_vf_with_tiles_2elem
    export blocking_vf_with_tiles_2elem_tiles

    include("./therm2D.jl")
    export sigma
    export set_bc_part!
    export tempsolver
    


end