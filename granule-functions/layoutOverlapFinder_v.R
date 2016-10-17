# Import the relevant packages (All for compiling the C++ code inline)
library(Rcpp)
library(RcppArmadillo)
library(inline)

# We need to include these namespaces in the C++ code 
includes <- '
using namespace Rcpp;
using namespace arma;
'

# This is the main C++ function 
# We cast 'm' as an Armadillo matrix 'm1' and compute the number of rows 'numRows'
# We cast 'x' as a row vector 'x1'
# We then loop through the rows of the matrix 
# As soon as we find a matching row (anyEqual = TRUE), we stop and return TRUE
# If no matching row is found, then anyEqual = FALSE and we return FALSE
# Note: Within the for loop, we do an elementwise comparison of a row of m1 to x1
# If the row is equal to x1, then the sum of the elementwise comparision should equal the number of elements of x1
src <- '
mat m1 = as<mat>(m); 
int numRows = m1.n_rows;
rowvec x1 = as<rowvec>(x);
bool anyEqual = FALSE;
for (int i = 0; i < numRows & !anyEqual; i++){
anyEqual = (sum(m1.row(i) == x1) == x1.size());
}
return(wrap(anyEqual));
'

# Here, we compile the function above
# Do this once (in a given R session) and use it as many times as desired
layoutOverlapFinder.v <- cxxfunction(signature(m="numeric", x="numeric"), src, plugin='RcppArmadillo', includes)





# We need to include these namespaces in the C++ code 
includes <- '
using namespace Rcpp;
using namespace arma;
'

# This is the main C++ function 
# We cast 'm' as an Armadillo matrix 'm1' and compute the number of rows 'numRows'
# We cast 'x' as a row vector 'x1'
# We then loop through the rows of the matrix 
# As soon as we find a matching row (anyEqual = TRUE), we stop and return TRUE
# If no matching row is found, then anyEqual = FALSE and we return FALSE
# Note: Within the for loop, we do an elementwise comparison of a row of m1 to x1
# If the row is equal to x1, then the sum of the elementwise comparision should equal the number of elements of x1
src <- '
mat m1 = as<mat>(m); 
int numRows = m1.n_rows;
mat x1 = as<mat>(x);
vec anyEqual = zeros<vec>(x1.n_rows);
for (int j = 0; j < x1.n_rows; j++){
for (int i = 0; i < numRows & !anyEqual(j); i++){
anyEqual(j) = (sum(m1.row(i) == x1.row(j)) == x1.n_cols);
}
}
return(wrap(anyEqual));
'

# Here, we compile the function above
# Do this once (in a given R session) and use it as many times as desired
layoutOverlapFinder.m <- cxxfunction(signature(m="numeric", x="numeric"), src, plugin='RcppArmadillo', includes)


