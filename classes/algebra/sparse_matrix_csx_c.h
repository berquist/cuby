//##############################################################################
// Sparse matrix data structure
//##############################################################################

// Prototypes for sparse_matrix_csx_c.c

extern SparseMatrixCSXData *SparseMatrixCSXData_new();
extern void SparseMatrixCSXData_free(SparseMatrixCSXData* md);
extern void SparseMatrixCSXData_alloc(SparseMatrixCSXData *ptr, int row_mode, int m, int n, int nnz, int recs);

// Helpers
extern void sparse_matrix_raise_index_error(SparseMatrixCSXData *a, int i, int j);
