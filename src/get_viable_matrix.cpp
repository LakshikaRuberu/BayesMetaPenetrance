#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.get_viable_matrix)]]

NumericMatrix get_viable_matrix(NumericMatrix input_matrix){
   
   int num_of_rows=input_matrix.nrow();
   int num_of_cols=input_matrix.ncol();     
  // Create a column matrix with all zeros
  //Its meant to hold 1 for rows satisfying the underlying condition
  // of x > x++ 
  
 NumericMatrix get_viable_matrix(
          NumericMatrix input_matrix){
   
   int num_of_rows=input_matrix.nrow();
   int num_of_cols=input_matrix.ncol();     
  // Create a column matrix with all zeros
  //Its meant to hold 1 for rows satisfying the underlying condition
  // of x > x++ 
  
  NumericMatrix viable_rows( num_of_rows  ,1);
  
  
  int xsize = viable_rows.nrow() * viable_rows.ncol(); // fill with 1
    for (int i = 0; i < xsize; i++) {
    viable_rows[i] = 1;
     }

    
  // Going Through All Rows 
  for (int row_index=0; row_index < num_of_rows; row_index++) 
  {
    // Getting Current Row 
    NumericMatrix::Row current_row = input_matrix( row_index , _ ); 
    
    for (int col_index=0; col_index < num_of_cols-1; col_index++) 
    {
      //checking condition 
      if(current_row[col_index]>=current_row[col_index+1])
      {
        viable_rows[row_index] = 0;      
        break;
      }// end if 
      
    }//end colums
  
  
  }//end row

  return viable_rows;

}

