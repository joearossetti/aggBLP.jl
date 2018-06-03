## example of loading data and working with it

## reading in data matrices
X1 = readcsv("inst/blp_test_data/X1_mat.csv", header=true)
X2 = readcsv("inst/blp_test_data/X2_mat.csv", header=true)
Z = readcsv("inst/blp_test_data/Z_mat.csv", header=true)
y = readcsv("inst/blp_test_data/y_mat.csv", header=true)
ids = readcsv("inst/blp_test_data/id_mat.csv", Int8, header=true)
true_delta = readcsv("inst/blp_test_data/true_delta_mat.csv", header=true)
true_sigma = readcsv("inst/blp_test_data/true_sigma_mat.csv", header=true)
true_alpha_beta = readcsv("inst/blp_test_data/true_alpha_beta_vec.csv", header=true)

true_sigma_parmtrs = diag(true_sigma[1])

true_sigma_mat = diagm(true_sigma_parmtrs)

using aggBLP

my_nsfp_state = aggBLP.BlpNsfpState(true_delta[1],
true_alpha_beta[1],
true_sigma_parmtrs,
true_sigma_mat,
nothing,
nothing,
nothing,
nothing)

dim_full_nln_mat = 4
nln_mat = my_nsfp_state.nln_parmtrs_mat
nln_parmtrs_vec = my_nsfp_state.nln_parmtrs
full_nln_parmtrs_vec = vec(nln_mat)

function short_to_long!(nsfp_state)
    nsfp_state.



## should be set at construction
num_nln_parmtrs = size(nln_mat, 1)
num_draws = 1000

using Sobol
my_sobol_sequence = SobolSeq(num_nln_parmtrs)
skip(my_sobol_sequence, num_draws)

my_sobol_sample = hcat([next(my_sobol_sequence) for j = 1:num_draws]...)
## end should be set at construction

chol(nln_mat)

## function to take a matrix and compute the determinants of its
## principal minors

function detprncplminor(X)
    @assert size(X)[1] == size(X)[2]

    dimX = size(X)[1]

    dets = Array{Float64}(dimX)

    for j = 1:dimX
        dets[j] = det(X[1:j,1:j])
    end

    return dets
end


## function to compute the adjugate and use jacobi's rule to get the
## derivatives of the determinants
## vectorize these matrices to get
## gradient of determinants of principal minors
## note that the first principal minor depends only on the first element so
## add 0s bordering its adjugate to get a nxn matrix

function grad_detprncplminor(X)
    @assert size(X)[1] == size(X)[2]

    dimX = size(X)[1]

    grad_dets = zeros(Array{Float64}(dimX, dimX^2))

    for j in 1:dimX
        adj_minor = inv(X[1:j,1:j]) * det(X[1:j,1:j])
        for h in 1:j
            grad_dets[j,h] = vec(adj_minor)[h]
        end
    end

    return grad_dets
end

## constraint function
function nln_mat_constraint(long_nln_parmtrs, dim_sigma_mat)
    return detprncplminor(reshape(long_nln_parmtrs, dim_sigma_mat, dim_sigma_mat))
end

## gradient of constraint
function grad_nln_mat_constraint(long_nln_parmtrs, dim_sigma_mat)
    return grad_detprncplminor(reshape(long_nln_parmtrs, dim_sigma_mat, dim_sigma_mat))
end
