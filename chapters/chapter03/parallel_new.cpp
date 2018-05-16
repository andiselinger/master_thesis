MatMatMultSymbolic_MPIAIJ_MPIAIJ_seq(Mat A, Mat B, Mat *C)
{
  // get B_non_local by taking rows of B (= non-zero columns of local A) from other processors
  MatGetBrowsOfAoCols_MPIAIJ(A, B, B_nonloc);
  // Symbolic calculation a_diag_b_diag = A_loc_diag * B_loc_diag
  MatMatMultSymbolic_SeqAIJ_SeqAIJ(A_loc->diag, B_loc->diag, &a_diag_b_diag);
  // Symbolic calculation of the A_loc_diag * B_loc_off
  for (i : all local diagonal rows of A) {
    for (j : all non zero off-diagonal columns of i) {
      // Correct the indices
      for(k : all non-zero columns in row j of B) {
        columns_b_row_loc_off[k] = B->offdiag->global_index[column_index[k]];
      }
      // Add non-zero local offdiagonal columns of B into the sorted linked list
      addColumnsToSortedList(columns_b_row_loc_off, num_columns_in_b_row, list);
    }
    updateIArray(num_cols, ci);
    // Copy the column indices from the sorted list to memspace
    copyListToMemSpace(list, mem_space);
    // Initialize the list with zeros
    reinitialize(list);
  }
  // Copy from linked list to j-array
  PetscFreeSpaceContiguous(&mem_space, a_loc_diag_b_loc_off);
  // Symbolic calc of A_local_off * B_nonloc
  MatMatMultSymbolic_SeqAIJ_SeqAIJ(a_loc_off, b_nonloc, &a_loc_off_b_nonloc);
  for (i : all local rows of A) {
    // Correct indices from A_loc_diag * B_loc_diag
    addRowStart(a_loc_diag_b_loc_diag);
    // Merge j-arrays of A_loc_diag * B_loc_diag  and A_loc_diag * B_loc_off and A_loc_off * B_nonloc
    Merge3SortedArrays(a_loc_diag_b_loc_diag,
                       a_loc_diag_b_loc_off,
                       a_loc_off_b_nonloc,
                       cj);
    // Add information about the number of non-zero columns in the current row to the array ci
    // Set the arrays dnz and onz which are needed for preallocation.
    // They give information about the number of zero in each row 
    // for diagonal and off-diagonal entries
    setDnzOnz(mem_space, dnz, onz);
    updateArray(num_columns, ci);
  }
  // Create and assemble parallel (multi-processor) matrix Cmpi
  MatCreate(&Cmpi);
  MatPreallocate(dnz,onz);
  for (i : all local rows of A) {
    cnz = getNumColumnsinRow(ci, i);
    MatSetValues(Cmpi, global_row_number, cnz, cj, ca);
  }
  MatAssembly(Cmpi);
}
