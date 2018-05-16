MatMatMultSymbolic_MPIAIJ_MPIAIJ_nonscalable(Mat A, Mat B, Mat *C)
{
  // get the non local part B_non_local by taking rows of B (= non-zero columns of local A) from other processors
  MatGetBrowsOfAoCols(A, B, B_non_local);
  // Compute symbolic A_loc * B = A_loc_diag*B_loc + A_loc_off*B_nonloc
  for (i : all local rows of A) {
    // diagonal portion of A
    for (j : all non-zero diagonal columns of i) {
      // Get non-zero columns of the current row of B
      columns_b_row = getRow(B_loc, j);
      // add non-zero columns of columns_b_row into the sorted linked list
      addColumnsToSortedList(columns_b_row,
                             num_columns_in_b_row,
                             list);
    }
    // off-diagonal portion of A
    for (j : all non-zero off diagonal columns of A) {
      // Get non-zero columns of the current row of B
      columns_b_row = getRow(B_nonloc, j);
      // add non-zero columns of columns_b_row into the sorted linked list
      addColumnsToSortedList(columns_b_row,
                             num_columns_in_b_row,
                             list);
    }
    // Add information about the number of non-zero columns in the current row to the array ci
    updateIArray(num_columns, ci);
    // Copy the column indices from the sorted list to the memory space
    copyListToMemSpace(list, mem_space);
    // Set the arrays dnz and onz which are needed for preallocation.
    // They give information about the number of zero in each row 
    // for diagonal and off-diagonal entries
    setDnzOnz(mem_space, dnz, onz);
    // Reinitialze the list so it is ready for a new row
    reinitializeList(list);
  }
  // Copy the results to the result j-array
  copyMemSpaceToJArray(mem_space, cj);
  // create and assemble symbolic parallel (multi-processor) matrix Cmpi
  createMatrix(Cmpi);
  // Preallocate the matrix with the arrays dnz and onz
  ierr = MatPreallocate(dnz,onz);
  for (i : all local rows of A) {
    cnz = getNumColumnsInRow(ci,i);
    MatSetValues(Cmpi,global_row_number,cnz,cj,ca);
  }
  MatAssembly(Cmpi);
}
