module aggBLP

mutable struct BlpNsfpState <: Any
    deltas
    ln_parmtrs
    nln_parmtrs
    gmm_obj_val
    cov_mat
end


end # module
