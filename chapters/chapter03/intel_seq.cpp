MatMatMult_SeqAIJ_SeqAIJ_combined(Mat A, Mat B, Mat *C)
{
  // Determine an upper bound on memory  requirements for C
  // Iterate over all rows of A (all rows of C)
  for (i : rows of A)
  {
    max_num_nonzeros = calcMaxNumNonZerosRow(B,i);
    c_maxmem += max_num_nonzeros;
  }
  // Allocate memory for the matrix C
  malloc(C);
  // Iterate over all rows of A (= over all rows of C)
  for (i : rows of A)
  {
    // Do the symbolic and numeric calculations
    // Iterate over all non zero values in the current row of A
    for (a_col : non-zero columns of i)
    {
      // Iterate over all non zeros in the row of B that
      // corresponds to the current column a_col
      for (b_col : non-zero columns of row in B)
      {
#IF APPEND_FLAG
        // If the flag of the current column is unset ...
        if (current_column_flag == UNSET)
        {
          // ... append the current column index to the
          // list of indices in this row of C.
          appendIndexToArray(current_column_index, c_row_j);
          // ... and set the corresponding flag
          setFlag(current_column_flag);
        }
#ELSE    // Look for index in c_row_j, if not present add it   
        insertInArray(current_column_index, c_row_j);
#ENDIF        
        // Add the product of:
        // - the value at the current column of A, and
        // - the value of the current column of B
        // to the dense array, at the current column of B index
        c_row_val_dense[b_col] += a_row_val[a_col] * b_row_val[b_col];
      }
    }
    // add information about the number of zeros
    // in the current row to the array ci
    updateIArray(ci);
    // Sort the array with the column indices of
    // the current row of C.
    sort(c_row_j);
    // Copy numeric values of the current row of C
    // from the dense array to the sparse array ca.
    copyDenseRowToSparse(c_row_val_dense, ca);
    // Unset the dense array with index flags of
    // the current row of C
    unsetIdxFlags(c_row_idx_flags);
    // Clear the numerical values of
    // the current row of C
    clearValues(c_row_val_dense);
  }   // End of the iteration over all rows of A
  // Create the new matrix from the arrays ca, ci, cj
  createSeqMatrix(C_seq);
}
