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

true_sigma_short = vec(true_sigma[1])[find(true_sigma[1].!=0.0)]

#dim_sigma = 4

using aggBLP
using Sobol

#indicator_mat = [1 0 0 0; 0 1 0 0; 0 0 1 0 ; 0 0 0 1]

#num_draws = 100

#init_short_pars = true_sigma_short

function update_nln_par_mat!(unobsv_hetero, short_par_vec)
      unobsv_hetero.nln_parmtrs_long[unobsv_hetero.which_nln_parmtrs] = short_par_vec
      return nothing
end

function update_nui_mat!(unobsv_hetero)
      print(unobsv_hetero.nln_parmtrs_mat)
      chol_decomp = chol(Hermitian(unobsv_hetero.nln_parmtrs_mat))
      unobsv_hetero.nui_mat = chol_decomp * unobsv_hetero.rand_mat
      return nothing
end


function create_uhobj(indicator_mat, init_short_pars, num_draws,dim_sigma)
      ## make the which_nln_parmtrs
      which_nln_parmtrs = find(indicator_mat .== 1)

      ## how many nln parameters?
      num_nln_parmtrs = length(init_short_pars)

      ## create rand_mat
      sobol_sequence = Sobol.SobolSeq(dim_sigma)
      Sobol.skip(sobol_sequence, num_draws)

      sobol_sample = hcat([Sobol.next(sobol_sequence) for j = 1:num_draws]...)

      nln_parmtrs_long = Array{Float64}(dim_sigma^2)

      nln_parmtrs_mat = reshape(nln_parmtrs_long, dim_sigma, dim_sigma)

      #nln_parmtrs_long[my_which_nln_pars] = init_short_pars

      ## initialize UnobsvHetero
      uhobj = aggBLP.UnobsvHetero(dim_sigma,
       num_draws,
       num_nln_parmtrs,
       which_nln_parmtrs,
       nln_parmtrs_long, # nln_parmtrs_long
       nln_parmtrs_mat, # nln_parmtrs_mat
       sobol_sample, # rand_mat
       Array{Float64}(dim_sigma, num_draws)) # nui_mat

       update_nln_par_mat!(uhobj, init_short_pars)
       update_nui_mat!(uhobj)

       return uhobj
end

# uhobj.nln_parmtrs_long[uhobj.which_nln_parmtrs] = init_short_pars
# chol_decomp = chol(uhobj.nln_parmtrs_mat)
# uhobj.nui_mat = chol_decomp * uhobj.rand_mat

my_uhobj = create_uhobj([1 0 0 0; 0 1 0 0; 0 0 1 0 ; 0 0 0 1], true_sigma_short, 100, 4)
