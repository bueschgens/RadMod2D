module RadMod2D

    # julia packages
    using DelimitedFiles
    using StaticArrays
    using LinearAlgebra
    using GLMakie
    # using CairoMakie # for image export
    # CairoMakie.activate!(type = "svg")
    
    include("./struct2D.jl")
    export Point2D, Index2D
    export norm, get_length

    include("./mesh2D.jl")
    export AbstractModel, MutableModel, ConstModel
    export create_empty_model
    export add!, offset_model!
    export ConstModel
    export element_analysis

    include("./geom2D.jl")
    export edge, rectangle, circle, circle_open, cosinus

    include("./raycast2D.jl")
    export create_randomly_occupied_tiles
    export get_max_steps
    export tilewalk_with_check_for_occ_tiles
    export tilewalk_with_return, tilewalk_with_return!

    include("./view2D.jl")
    export existing_vf!
    export blocking_vf_brute_force!
    export blocking_vf_with_tiles!
    export blocking_vf_with_tiles_simplified!
    export calculating_vf!
    export get_area_of_part
    export get_tile_dimensions


    include("./view2Dchecks.jl")
    # comapre elements checks
    export blocking_vf_with_tiles_2elem
    export blocking_vf_with_tiles_2elem_tiles
    export check_tile_occupation
    export check_tile_occupation_parts
    export tile_occ_analysis

    include("./view2Dpostproc.jl")
    export compact_vfmat_to_parts
    export save_vfmat_part_to_csv

    include("./therm2D.jl")
    export sigma
    export set_bc_part!
    export tempsolver

    include("./plot2D.jl")

    export new_figure
    export plot_model, plot_model_elements, plot_empty_tiles, plot_occupied_tiles, 
        plot_points_and_connection, plot_existing, plot_model_shadow_1to1, 
        plot_model_shadow_1toAll, plot_model_with_value, setup_axis!

end