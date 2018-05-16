MatMatMultNumeric_MPIAIJ_MPIAIJ(Mat A,Mat P,Mat C)
{
  // get B_non_local by taking rows of B (= non-zero columns of local A) from 
  // other processors
  MatGetBrowsOfAoCols_MPIAIJ(A, B, B_nonlocal);
  // Iterate over all local rows of C
  for (i : all local rows of C) {
    // Calculate A_loc_diag * B_loc
    for (j : all non-zero diagonal columns of A) {
      for (k : all non-zero columns of row i in C) {
        if (k is non-zero in row j of B) {
          c_sparse[k] += a[j]*b[k];
        }
      }
    }
    // Calculate A_loc_off * B_nonloc
    for (j : all non-zero diagonal columns of A) {
      for (k : all non-zero columns of row i in  C) {
        if (k is non-zero in row j of B) {
          c_sparse[k] += a[j]*b[k];
        }
      }
    }

    // Copy the results from the sparse row to the off-diagonal and diagonal
    // portions of the result matrix C
    // 1st off-diagonal part of C
    for (k : non-zero off-diagonal columns of i in C that precede the diagonal) 
    {
      c_off[correctIndex(k)] = c_sparse[k];
      c_sparse[k] = 0.0;
    }

    // diagonal part of C
    for (k : all non-zero diagonal columns of i in C) {
      c_diag[correctIndex(k)] = c_sparse[k];
      c_sparse[k] = 0.0;
    }

    // 2nd off-diagoanl part of C
    for (k : non-zero off-diagonal columns of i in C that succeed the diagonal) 
    {
      c_off[correctIndex(k)] = c_sparse[k];
      c_sparse[k] = 0.0;
    }
  }
  // Finalize data structure of C
  ierr = MatAssembly(C);
}
