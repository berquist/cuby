#ifndef ALGEBRA_C_H
#define ALGEBRA_C_H

//==============================================================================
// Compatibility fix for ruby version < 1.8.5
//==============================================================================
#ifndef RARRAY_PTR
        #define RARRAY_PTR(v) RARRAY(v)->ptr
#endif
#ifndef RARRAY_LEN
        #define RARRAY_LEN(v) RARRAY(v)->len
#endif
#ifndef RSTRING_PTR
        #define RSTRING_PTR(v) RSTRING(v)->ptr
#endif
#ifndef RSTRING_LEN
        #define RSTRING_LEN(v) RSTRING(v)->len
#endif


//==============================================================================
// Ruby classes
//==============================================================================
extern VALUE rb_mAlgebra;
extern VALUE rb_eDimensionError;
extern VALUE rb_eSparseWriteError;
extern VALUE rb_eSparseTypeError;
extern VALUE rb_cMatrix;
extern VALUE rb_cVector;
extern VALUE rb_cSparseMatrixCSX;

//==============================================================================
// Type used for matrix and vector elements
//==============================================================================
typedef double dbl;

//==============================================================================
// Matrix data structure
//==============================================================================
// m by n matrix: m rows, n columns
// indexed as Aij: i is column, j is row
// internally stored in fortran format
typedef struct matrixdata {
	int m;
	int n;
	dbl *data;
} MatrixData;

//==============================================================================
// Vector data structure
//==============================================================================
// vector of size n
typedef struct vectordata {
	int size; // n
	dbl *data;
} VectorData;

//==============================================================================
// Sparse matrix data structure
//==============================================================================
// sparse matrix in CSR format
typedef struct sparsematrixcsxdata {
	int row_mode;	// CSX mode, 1 for row, 0 for col
	int m;		// number of rows
	int n; 		// number of columns
	int nnz;	// number of nonzero elements
	dbl *data_values;
	int *data_records;
	int *data_indices;
} SparseMatrixCSXData;

//==============================================================================
// Prototypes: Init methods
//==============================================================================
void Init_matrix_c ();
void Init_vector_c ();
void Init_sparse_matrix_csx_c ();

//==============================================================================
// Prototypes: Vector data struct
//==============================================================================
VectorData *VectorData_new();
void VectorData_free(VectorData* md);
void VectorData_alloc(VectorData *ptr, int n);
void VectorData_alloc_zero(VectorData *ptr, int n);

VALUE rb_cVector_new(VALUE class);

// Vector type check:
void vector_raise_if_not_vector(VALUE val);

//==============================================================================
// Prototypes: Matrix data struct
//==============================================================================
MatrixData *MatrixData_new();
void MatrixData_free(MatrixData* md);
void MatrixData_alloc(MatrixData *ptr, int m, int n);
void MatrixData_alloc_zero(MatrixData *ptr, int m, int n);

VALUE rb_cMatrix_new(VALUE class);

//==============================================================================
// Prototypes: SparseMatrixCSX data struct
//==============================================================================
SparseMatrixCSXData *SparseMatrixCSXData_new();
void SparseMatrixCSXData_free(SparseMatrixCSXData* md);
void SparseMatrixCSXData_alloc(SparseMatrixCSXData *ptr, int row_mode, int m, int n, int nnz, int rec);

VALUE rb_cSparseMatrixCSX_new(VALUE class);

#endif
