library(Rcpp)
cppFunction('NumericVector matrixIndex(NumericMatrix m1, NumericMatrix m2){

            int m1Rows = m1.nrow();
            int m2Rows = m2.nrow();
            NumericVector out(m2Rows);  
            
            for (int i = 0; i < m1Rows; i++){
            for (int j = 0; j < m2Rows; j++){
            
            if(m1(i, 0) == m2(j, 0) && m1(i, 1) == m2(j, 1)){
            out[j] = (j+1);
            //out.push_back(j + 1);
            }
            }
            }
            
            return out;
            
            }')