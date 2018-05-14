## example of loading data and working with it

## reading in data matrices

X1 = readcsv("inst/blp_test_data/X1_mat.csv", header=true)
X2 = readcsv("inst/blp_test_data/X2_mat.csv", header=true)
Z = readcsv("inst/blp_test_data/Z_mat.csv", header=true)
y = readcsv("inst/blp_test_data/y_mat.csv", header=true)
ids = readcsv("inst/blp_test_data/id_mat.csv", Int8, header=true)
