module aggBLP

mutable struct UnobsvHetero <: Any
    dim_sigma
    num_draws
    num_nln_parmtrs
    which_nln_parmtrs
    nln_parmtrs_long
    nln_parmtrs_mat
    rand_mat
    nui_mat
end

mutable struct BlpNsfpState <: Any
    deltas
    ln_parmtrs
    unobsv_hetero
    pred_shares
end

mutable struct BlpNsfpObj <: Any
    nsfp_state::BlpNsfpState
    #useful_quantities e.g. pre-computing projection matrices and stuff
end

## see https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1
function (Q::BlpNsfpObj)(theta)
    ## in here do stuff to compute the gmm_obj function
    return theta
end

end # module
