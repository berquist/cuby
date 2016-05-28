#include <ruby.h>
#include <math.h>

#include "algebra_c.h"

#ifdef LAPACK_FOUND
#include "clapack_proto.h"
#endif


// Classes defined here
VALUE rb_cMatrix;

MatrixData *MatrixData_new() {
	MatrixData *ptr = (MatrixData*) malloc(sizeof(MatrixData));
	ptr->m = 0;
	ptr->n = 0;
	ptr->data = NULL;
	return ptr;
}

void MatrixData_free(MatrixData* md) {
	if (md->data != NULL) free(md->data);
	free(md);
}

void MatrixData_alloc(MatrixData *ptr, int m, int n) {
	if (m <= 0) rb_raise(rb_eTypeError, "Matrix dimension m must be larger than 0");
	if (n <= 0) rb_raise(rb_eTypeError, "Matrix dimension n must be larger than 0");
	ptr->m = m;
	ptr->n = n;
	ptr->data = (dbl*) malloc(m * n * sizeof(dbl));
}

void MatrixData_alloc_zero(MatrixData *ptr, int m, int n) {

	MatrixData_alloc(ptr, m, n);
	int size = m * n;
	int i;
	for(i = 0; i < size; i++) ptr->data[i] = 0.0;
}

MatrixData *MatrixData_transpose(MatrixData *ptr) {
	MatrixData *result = MatrixData_new();
	MatrixData_alloc(result, ptr->n, ptr->m);
	int i, j;
	for (i = 0; i < result->m; i++) {
		for(j = 0; j < result->n; j++) {
			result->data[i + j*result->m] = ptr->data[j + i*ptr->m];
		}
	}
	return result;
}

void MatrixData_fill(MatrixData *ptr, dbl fill){
	int size = ptr->m * ptr->n;
	int i;
	for(i = 0; i < size; i++) ptr->data[i] = fill;
}

//##############################################################################
// Helpers
//##############################################################################

void matrix_raise_if_not_same_size(MatrixData *a, MatrixData *b) {
	if (a->n != b->n || a->m != b->m) {
		rb_raise(rb_eDimensionError, "Matrices must be of same size");
	}
}

void matrix_raise_if_not_square(MatrixData *a) {
	if (a->n != a->m ) {
		rb_raise(rb_eDimensionError, "Matrix must be square");
	}
}

void matrix_raise_index_error(MatrixData *a, int i, int j) {
	if (i >= a->m || i < 0) {
		rb_raise(rb_eIndexError, "Row index (i) out of bounds");
	}
	if (j >= a->n || j < 0) {
		rb_raise(rb_eIndexError, "Column index (j) out of bounds");
	}
}

void matrix_raise_if_not_matrix(VALUE val) {
	if (rb_obj_is_kind_of(val, rb_cMatrix) != Qtrue) {
		rb_raise(rb_eTypeError, "Argument is not a matrix");
	}
}

//##############################################################################
// Matrix class methods
//##############################################################################

//==============================================================================
// Allocation function
//==============================================================================

VALUE rb_cMatrix_new(VALUE class) {
	MatrixData *ptr = MatrixData_new();
	VALUE data = Data_Wrap_Struct(class, 0, MatrixData_free, ptr);
	return data;
}

//==============================================================================
// Constructors
//==============================================================================

VALUE rb_cMatrix_initialize(VALUE self, VALUE colsize, VALUE rowsize) {
	MatrixData *ptr;
	int m = NUM2INT(colsize);
	int n = NUM2INT(rowsize);
	Data_Get_Struct(self,MatrixData, ptr);
	MatrixData_alloc_zero(ptr, m, n);
	return Qnil;
}

VALUE rb_cMatrix_filled(VALUE class, VALUE colsize, VALUE rowsize, VALUE fill) {
	int m = NUM2INT(colsize);
	int n = NUM2INT(rowsize);
	dbl fillnum = NUM2DBL(fill);
	VALUE result = rb_cMatrix_new(class);
	MatrixData *ptr;
	Data_Get_Struct(result,MatrixData, ptr);
	MatrixData_alloc(ptr, m, n);
	MatrixData_fill(ptr, fillnum);
	return result;
}

VALUE rb_cMatrix_zero(int argc, VALUE *argv, VALUE class) {
	// Parse arguments
	VALUE v_m, v_n;
	rb_scan_args(argc, argv, "11", &v_m, &v_n);
	int m = NUM2INT(v_m);
	int n;
	if (v_n == Qnil) n = m; else n = NUM2INT(v_n);
	VALUE result = rb_cMatrix_new(class);
	MatrixData *ptr;
	Data_Get_Struct(result,MatrixData, ptr);
	MatrixData_alloc_zero(ptr, m, n);
	return result;
}

VALUE rb_cMatrix_from_array(VALUE class, VALUE array) {
	VALUE result = rb_cMatrix_new(class);
	MatrixData *ptr;
	Data_Get_Struct(result,MatrixData, ptr);
	Check_Type(array, T_ARRAY); // Check if array is Array
	int csize = RARRAY_LEN(array);
	if (csize == 0) rb_raise(rb_eTypeError, "Array can not be empty.");
	VALUE* array_ptr = RARRAY_PTR(array);
	VALUE row = array_ptr[0];
	VALUE* row_ptr;
	Check_Type(row, T_ARRAY); // Check if first row is Array
	int rsize = RARRAY_LEN(row);
	if (rsize == 0) rb_raise(rb_eTypeError, "Array can not be empty.");
	MatrixData_alloc(ptr, csize, rsize);

	int i, j;
	for (j = 0; j < csize; j++) {
		row = array_ptr[j];
		// Check if row is Array
		Check_Type(row, T_ARRAY);
		// Check row for size
		if (RARRAY_LEN(row) != rsize) rb_raise(rb_eTypeError, "Rows must have the same size!");
		for(i = 0; i < rsize; i++) {
			row_ptr = RARRAY_PTR(row);
			ptr->data[j + i*csize] = NUM2DBL(row_ptr[i]);
		}
	}
	return result;
}

// diagonal matrix with specified number on diagonal
VALUE rb_cMatrix_diagonal_1(VALUE class, VALUE size, VALUE number) {
	int n = NUM2INT(size);
	dbl fill = NUM2DBL(number);
	// Alocate result
	VALUE result = rb_cMatrix_new(class);
	MatrixData *ptr_r;
	Data_Get_Struct(result,MatrixData, ptr_r);
	MatrixData_alloc_zero(ptr_r, n, n);
	int i;
	for (i = 0; i < n; i++) {
		ptr_r->data[i + i* ptr_r->m] = fill;
	}
	return result;
}

// diagonal matrix with values from an array on the diagonal
VALUE rb_cMatrix_diagonal_2(VALUE class, VALUE array) {
	int n = RARRAY_LEN(array);
	// Alocate result
	VALUE result = rb_cMatrix_new(class);
	MatrixData *ptr_r;
	Data_Get_Struct(result,MatrixData, ptr_r);
	MatrixData_alloc_zero(ptr_r, n, n);
	VALUE* array_ptr;
	int i;
	for (i = 0; i < n; i++) {
		array_ptr = RARRAY_PTR(array);
		ptr_r->data[i + i* ptr_r->m] = NUM2DBL(array_ptr[i]);
	}
	return result;

}

//==============================================================================
// Matrix information
//==============================================================================

VALUE rb_cMatrix_size_m(VALUE self) {
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	return INT2NUM(ptr->m);
}

VALUE rb_cMatrix_size_n(VALUE self) {
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	return INT2NUM(ptr->n);
}

//==============================================================================
// Elemet access
//==============================================================================

VALUE rb_cMatrix_get_ij(VALUE self, VALUE i, VALUE j) {
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	int ii = NUM2INT(i);
	int jj = NUM2INT(j);
	// Out of bounds check
	matrix_raise_index_error(ptr, ii, jj);
	return rb_float_new(ptr->data[ii + jj* ptr->m]);
}

VALUE rb_cMatrix_set_ij(VALUE self, VALUE i, VALUE j, VALUE val) {
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	int ii = NUM2INT(i);
	int jj = NUM2INT(j);
	// Out of bounds check
	matrix_raise_index_error(ptr, ii, jj);
	ptr->data[ii + jj* ptr->m] = NUM2DBL(val);
	return val;
}

VALUE rb_cMatrix_submatrix(VALUE self, VALUE i, VALUE j, VALUE m, VALUE n) {
	// Get parameters
	int ii = NUM2INT(i);
	int jj = NUM2INT(j);
	int mm = NUM2INT(m);
	int nn = NUM2INT(n);
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Check parameters
	matrix_raise_index_error(ptr, ii, jj);
	matrix_raise_index_error(ptr, ii+mm-1, jj+nn-1);
	/// Size 0 check
	// Allocate result
	VALUE result = rb_cMatrix_new(rb_class_of(self));
	MatrixData *ptr_r;
	Data_Get_Struct(result,MatrixData, ptr_r);
	MatrixData_alloc(ptr_r, mm, nn);
	// Copy data
	int a,b;
	for (a=0; a<mm; a++) {
		for (b=0; b<nn; b++) {
			ptr_r->data[a + b * ptr_r->m] = ptr->data[a+ii + (b+jj) * ptr->m];
		}
	}	
	// Return
	return result;
}

//==============================================================================
// Row & col access
//==============================================================================

VALUE rb_cMatrix_col_to_a(VALUE self, VALUE index) {
	int i, j;
	j = NUM2INT(index);
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	matrix_raise_index_error(ptr, 0, j);
	// New array
	VALUE array = rb_ary_new2(ptr->m);
	// Fill array
	for (i = 0; i < ptr->m; i++) {
		rb_ary_store(array, i, rb_float_new(ptr->data[i + j * ptr->m])); 
	}
	return array;
}

VALUE rb_cMatrix_row_to_a(VALUE self, VALUE index) {
	int i, j;
	i = NUM2INT(index);
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	matrix_raise_index_error(ptr, i, 0);
	// New array
	VALUE array = rb_ary_new2(ptr->n);
	// Fill array
	for (j = 0; j < ptr->n; j++) {
		rb_ary_store(array, j, rb_float_new(ptr->data[i + j * ptr->m])); 
	}
	return array;
}

VALUE rb_cMatrix_col_to_vector(VALUE self, VALUE index) {
	int i, j;
	j = NUM2INT(index);
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	matrix_raise_index_error(ptr, 0, j);
	// Allocate new vector
	VALUE result = rb_cVector_new(rb_cVector);
	VectorData *v_ptr;
	Data_Get_Struct(result,VectorData, v_ptr);
	VectorData_alloc(v_ptr, ptr->m);
	// Copy
	for (i = 0; i < ptr->m; i++) {
		v_ptr->data[i] = ptr->data[i + j * ptr->m];
	}
	return result;
}

VALUE rb_cMatrix_row_to_vector(VALUE self, VALUE index) {
	int i, j;
	i = NUM2INT(index);
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	matrix_raise_index_error(ptr, i, 0);
	// Allocate new vector
	VALUE result = rb_cVector_new(rb_cVector);
	VectorData *v_ptr;
	Data_Get_Struct(result,VectorData, v_ptr);
	VectorData_alloc(v_ptr, ptr->n);
	// Copy
	for (j = 0; j < ptr->n; j++) {
		v_ptr->data[j] = ptr->data[i + j * ptr->m];
	}
	return result;
}

VALUE rb_cMatrix_diagonal_to_vector(VALUE self) {
	int i, n;
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Get smaller dimension
	n = ptr->n;
	if (ptr->m < ptr->n) n = ptr->m;
	// Allocate new vector
	VALUE result = rb_cVector_new(rb_cVector);
	VectorData *v_ptr;
	Data_Get_Struct(result,VectorData, v_ptr);
	VectorData_alloc(v_ptr, n);
	// Copy
	for (i = 0; i < n; i++) {
		v_ptr->data[i] = ptr->data[i + i * ptr->m];
	}
	return result;
}


//==============================================================================
// Misc
//==============================================================================

VALUE rb_cMatrix_to_s(VALUE self) {
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	//printf("Size: %d x %d\n", ptr->n, ptr->m);
	int i, j;
	char buff[50];
	VALUE str = rb_str_buf_new(100);
	for (j=0; j < ptr->m; j++) {
		sprintf(buff,"[ ");
		rb_str_cat2(str,buff);
		for (i = 0; i < ptr->n; i++) {
			sprintf(buff,"%8.3f ",ptr->data[j + i*ptr->m]);
			rb_str_cat2(str,buff);
		}
		sprintf(buff,"]\n");
		rb_str_cat2(str,buff);
	}

	return str;
}

VALUE rb_cMatrix_to_a(VALUE self) {
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	VALUE row;
	VALUE array = rb_ary_new2(ptr->m);
	int i, j;
	for (i = 0; i < ptr->m; i++) {
		row = rb_ary_new2(ptr->n);
		for (j = 0; j < ptr->n; j++) {
			rb_ary_store(row, j, rb_float_new(ptr->data[i + j * ptr->m])); 
		}
		rb_ary_store(array, i, row);
	}
	return array;
}

VALUE rb_cMatrix_to_vector(VALUE self) {
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Check size - only column matrix permitted
	if (ptr->n != 1) rb_raise(rb_eDimensionError, "Matrix.to_vector works only with single column matrices");
	// Allocate new vector
	VALUE result = rb_cVector_new(rb_cVector);
	VectorData *v_ptr;
	Data_Get_Struct(result,VectorData, v_ptr);
	VectorData_alloc(v_ptr, ptr->m);
	// Copy
	memcpy(v_ptr->data, ptr->data, ptr->m * sizeof(dbl));
	return result;
}

VALUE rb_cMatrix_deep_copy(VALUE self){
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Alocate result
	VALUE result = rb_cMatrix_new(rb_class_of(self));
	MatrixData *ptr_r;
	Data_Get_Struct(result,MatrixData, ptr_r);
	MatrixData_alloc(ptr_r, ptr->m, ptr->n);
	// Memcpy
	memcpy(ptr_r->data, ptr->data, ptr->m *  ptr->n * sizeof(dbl));
	return result;
}

//==============================================================================
// Operators
//==============================================================================

VALUE rb_cMatrix_matrix_multiply(VALUE self, VALUE matrix) {
	matrix_raise_if_not_matrix(matrix);
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	MatrixData *ptr_b;
	Data_Get_Struct(matrix,MatrixData, ptr_b);
	// Check sizes
	if (ptr_a->n != ptr_b->m) rb_raise(rb_eDimensionError, "Wrong size of matrices to multiply");

	// Alocate result
	VALUE result;
	result = rb_cMatrix_new(rb_class_of(self));
	MatrixData *ptr_r;
	Data_Get_Struct(result,MatrixData, ptr_r);
	MatrixData_alloc(ptr_r, ptr_a->m, ptr_b->n);

	// Multiply

#ifdef BLAS_FOUND
	/* using blas */
	// Calculates C = a*A*B + b*C
	// dimensions: A is m*k, B is k*n, so C has to be m*n
	int m = ptr_a->m;
	int n = ptr_b->n;
	int k = ptr_b->m;
	double alpha = 1.0;
	double beta = 0.0;
	dgemm_("N","N", &m, &n, &k, &alpha, ptr_a->data, &m, ptr_b->data, &k, &beta, ptr_r->data, &m);
#else
	/* in C  (about 3x slover than blas) */
	int i,j,k;
	double sum;
	for(i=0; i<ptr_a->m; i++){
		for(j=0; j<ptr_b->n; j++) {
			sum = 0.0;
			for(k=0; k<ptr_a->n; k++){
				sum += ptr_a->data[i + k* ptr_a->m] * ptr_b->data[k + j* ptr_b->m];
			}
			ptr_r->data[i + j* ptr_r->m] = sum;
		}
	}
#endif
	// Return
	return result;
}

VALUE rb_cMatrix_matrix_multiply_yield_diagonal(VALUE self, VALUE matrix) {
	matrix_raise_if_not_matrix(matrix);
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	MatrixData *ptr_b;
	Data_Get_Struct(matrix,MatrixData, ptr_b);
	// Check sizes
	if (ptr_a->n != ptr_b->m) rb_raise(rb_eDimensionError, "Wrong size of matrices to multiply");

	// Get smaller dimension
	int n;
	n = ptr_a->m;
	if (ptr_b->n < ptr_a->m) n = ptr_b->n;

	// Alocate result Vector
	VALUE result = rb_cVector_new(rb_cVector);
	VectorData *ptr_r;
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, n);

	// Multiply
	int i,k;
	double sum;
	for(i=0; i<n; i++){
		sum = 0.0;
		for(k=0; k<ptr_a->n; k++){
			sum += ptr_a->data[i + k* ptr_a->m] * ptr_b->data[k + i * ptr_b->m];
		}
		ptr_r->data[i] = sum;
	}
	// Return
	return result;
}


VALUE rb_cMatrix_matrix_vector_multiply(VALUE self, VALUE vector) {
	vector_raise_if_not_vector(vector);
	// Get matrix data
	MatrixData *ptr_m;
	Data_Get_Struct(self,MatrixData, ptr_m);
	// Get vector data
	VectorData *ptr_v;
	Data_Get_Struct(vector,VectorData, ptr_v);
	// Check size
	if (ptr_m->n != ptr_v->size) rb_raise(rb_eDimensionError, "Cant multiply matrix with vector of wrong size");
	// Result variables
	VALUE result;
	VectorData *ptr_r;
	int i,j;
	double sum;
	// The result is Vector
	result = rb_cVector_new(rb_class_of(vector));
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, ptr_m->m);
	// Multiply
	for (i=0; i<ptr_m->m; i++) {
		sum = 0.0;
		for (j=0; j<ptr_m->n; j++) {
			sum += ptr_m->data[i + j* ptr_m->m] * ptr_v->data[j];
		}
		ptr_r->data[i] = sum;
	}
	return result;
}

VALUE rb_cMatrix_plus(VALUE self, VALUE matrix) {
	matrix_raise_if_not_matrix(matrix);
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	MatrixData *ptr_b;
	Data_Get_Struct(matrix,MatrixData, ptr_b);
	// Check dimensions
	matrix_raise_if_not_same_size(ptr_a, ptr_b);
	// Allocate result
	VALUE result = rb_cMatrix_new(rb_class_of(self));
	MatrixData *ptr_r;
	Data_Get_Struct(result,MatrixData, ptr_r);
	MatrixData_alloc(ptr_r, ptr_a->m, ptr_a->n);
	// Sum
	int size = ptr_a->m * ptr_a->n;
	int i;
	for (i = 0; i < size; i++) {
		ptr_r->data[i] = ptr_a->data[i] + ptr_b->data[i];
	}
	// Return
	return result;
}

VALUE rb_cMatrix_minus(VALUE self, VALUE matrix) {
	matrix_raise_if_not_matrix(matrix);
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	MatrixData *ptr_b;
	Data_Get_Struct(matrix,MatrixData, ptr_b);
	// Check dimensions
	matrix_raise_if_not_same_size(ptr_a, ptr_b);
	// Allocate result
	VALUE result = rb_cMatrix_new(rb_class_of(self));
	MatrixData *ptr_r;
	Data_Get_Struct(result,MatrixData, ptr_r);
	MatrixData_alloc(ptr_r, ptr_a->m, ptr_a->n);
	// Sum
	int size = ptr_a->m * ptr_a->n;
	int i;
	for (i = 0; i < size; i++) {
		ptr_r->data[i] = ptr_a->data[i] - ptr_b->data[i];
	}
	// Return
	return result;
}

VALUE rb_cMatrix_plus_inplace(VALUE self, VALUE matrix) {
	matrix_raise_if_not_matrix(matrix);
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	MatrixData *ptr_b;
	Data_Get_Struct(matrix,MatrixData, ptr_b);
	// Check dimensions
	matrix_raise_if_not_same_size(ptr_a, ptr_b);
	// Sum
	int size = ptr_a->m * ptr_a->n;
	int i;
	for (i = 0; i < size; i++) {
		ptr_a->data[i] = ptr_a->data[i] + ptr_b->data[i];
	}
	// Return
	return Qnil;
}

VALUE rb_cMatrix_minus_inplace(VALUE self, VALUE matrix) {
	matrix_raise_if_not_matrix(matrix);
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	MatrixData *ptr_b;
	Data_Get_Struct(matrix,MatrixData, ptr_b);
	// Check dimensions
	matrix_raise_if_not_same_size(ptr_a, ptr_b);
	// Sum
	int size = ptr_a->m * ptr_a->n;
	int i;
	for (i = 0; i < size; i++) {
		ptr_a->data[i] = ptr_a->data[i] - ptr_b->data[i];
	}
	// Return
	return Qnil;
}

VALUE rb_cMatrix_comutative_multiply(VALUE self, VALUE number) {
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	double n = NUM2DBL(number);
	// Allocate result
	VALUE result = rb_cMatrix_new(rb_class_of(self));
	MatrixData *ptr_r;
	Data_Get_Struct(result,MatrixData, ptr_r);
	MatrixData_alloc(ptr_r, ptr_a->m, ptr_a->n);
	// Traverse matrix
	int size = ptr_a->m * ptr_a->n;
	int i;
	for (i = 0; i < size; i++) {
		ptr_r->data[i] = ptr_a->data[i] * n;
	}
	// Return
	return result;
}

VALUE rb_cMatrix_multiply_inplace(VALUE self, VALUE number) {
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	double n = NUM2DBL(number);
	// Traverse matrix
	int size = ptr_a->m * ptr_a->n;
	int i;
	for (i = 0; i < size; i++) {
		ptr_a->data[i] = ptr_a->data[i] * n;
	}
	// Return
	return Qnil;
}

//==============================================================================
// Comparison
//==============================================================================

VALUE rb_cMatrix_equal(VALUE self, VALUE matrix) {
	if (rb_obj_is_kind_of(matrix, rb_cMatrix) != Qtrue) return Qfalse;
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	MatrixData *ptr_b;
	Data_Get_Struct(matrix,MatrixData, ptr_b);
	// Check dimensions
	matrix_raise_if_not_same_size(ptr_a, ptr_b);
	// Sum
	int size = ptr_a->m * ptr_a->n;
	int i;
	for (i = 0; i < size; i++) {
		if (ptr_a->data[i] != ptr_b->data[i]) return Qfalse;
	}
	return Qtrue;
}

VALUE rb_cMatrix_equal_eps(VALUE self, VALUE matrix) {
	if (rb_obj_is_kind_of(matrix, rb_cMatrix) != Qtrue) return Qfalse;
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	MatrixData *ptr_b;
	Data_Get_Struct(matrix,MatrixData, ptr_b);
	// Check dimensions
	matrix_raise_if_not_same_size(ptr_a, ptr_b);
	// Sum
	int size = ptr_a->m * ptr_a->n;
	int i;
	dbl epsilon = NUM2DBL(rb_cv_get(rb_cMatrix, "@@default_epsilon"));
	for (i = 0; i < size; i++) {
		if (fabs(ptr_a->data[i] - ptr_b->data[i]) > epsilon) return Qfalse;
	}
	return Qtrue;
}

//==============================================================================
// Matrix operations
//==============================================================================

VALUE rb_cMatrix_transpose(VALUE self) {
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	MatrixData *ptr_t = MatrixData_transpose(ptr);
	VALUE result = Data_Wrap_Struct(rb_class_of(self), 0, MatrixData_free, ptr_t);
	return result;
}

VALUE rb_cMatrix_transpose_inplace(VALUE self) {
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	MatrixData *ptr_t = MatrixData_transpose(ptr);
	ptr->m = ptr_t->m;
	ptr->n = ptr_t->n;
	free(ptr->data);
	ptr->data = ptr_t->data;
	free(ptr_t);
	return Qnil;
}

#ifdef LAPACK_FOUND
// Return array of 3 matrices: Eigenvectors, real, imag
VALUE rb_cMatrix_eigensystem(VALUE self) {
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Check if is square
	matrix_raise_if_not_square(ptr);
	
	// Allocate results and temporary 
	// Matrix a, contains input
	VALUE a = rb_cMatrix_deep_copy(self);
	MatrixData *ptr_a;
	Data_Get_Struct(a,MatrixData, ptr_a);

	// Matrix real, imag
	VALUE real = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_r;
	Data_Get_Struct(real,MatrixData, ptr_r);
	MatrixData_alloc(ptr_r, ptr->m, 1);

	VALUE imag = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_i;
	Data_Get_Struct(imag,MatrixData, ptr_i);
	MatrixData_alloc(ptr_i, ptr->m, 1);

	// Matrix vr, eigenvectors
	VALUE vr = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_vr;
	Data_Get_Struct(vr,MatrixData, ptr_vr);
	MatrixData_alloc(ptr_vr, ptr->m, ptr->n);
	
	int info[1];
	int n = ptr->m;

	int workspace_size = 1;
	dbl *workspace;
	workspace = (dbl*) malloc(workspace_size * sizeof(dbl));

	int one = 1;
	int minus_one = -1;

	// query
	dgeev_("N", "V", &n, ptr_a->data, &n, ptr_r->data, ptr_i->data, NULL, &one, ptr_vr->data, &n, workspace, &minus_one, info);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack query for workspace failed");
	workspace_size = (int)*workspace; // size of new workspace
	free(workspace);

	// calculation
	workspace = (dbl*) malloc(workspace_size * sizeof(dbl));
	dgeev_("N", "V", &n, ptr_a->data, &n, ptr_r->data, ptr_i->data, NULL, &one, ptr_vr->data, &n, workspace, &workspace_size, info);
	free(workspace);

	// Build array with results
	VALUE array = rb_ary_new2(3);
	rb_ary_store(array, 0, vr);
	rb_ary_store(array, 1, real);
	rb_ary_store(array, 2, imag);

	return array;
}
#endif

#ifdef LAPACK_FOUND
VALUE rb_cMatrix_lu(VALUE self) {
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Matrix a, contains input
	VALUE a = rb_cMatrix_deep_copy(self);
	MatrixData *ptr_a;
	Data_Get_Struct(a,MatrixData, ptr_a);
	// Calculation
	int info[1];
	int* piv_array;
	int min_mn = ptr_a->m; if (ptr_a->n < ptr_a->m) min_mn = ptr_a->n;
	//int max_mn = ptr_a->m; if (ptr_a->n > ptr_a->m) max_mn = ptr_a->n;

	piv_array = (int*) malloc(min_mn * sizeof(int));
	dgetrf_(&(ptr_a->m), &(ptr_a->n), ptr_a->data, &(ptr_a->m), piv_array, info);
	/// Check errors
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack LU decomposition failed");
	// Allocate output matrix P
	VALUE p = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_p;
	Data_Get_Struct(p,MatrixData, ptr_p);
	MatrixData_alloc_zero(ptr_p, ptr_a->m, ptr_a->m);
	// Allocate output matrix L
	VALUE l = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_l;
	Data_Get_Struct(l,MatrixData, ptr_l);
	MatrixData_alloc_zero(ptr_l, ptr->m, min_mn);
	// Allocate output matrix U
	VALUE u = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_u;
	Data_Get_Struct(u,MatrixData, ptr_u);
	MatrixData_alloc_zero(ptr_u, min_mn, ptr->n);

	// Copy L and U
	int i,j;
	for (i=0; i<ptr_a->m; i++) {
		for (j=0; j<ptr_a->n; j++) {
			if (j >= i) ptr_u->data[i + j * ptr_u->m] = ptr_a->data[i + j * ptr_a->m];
			if (j == i) ptr_l->data[i + j * ptr_l->m] = 1.0;
			if (j < i)  ptr_l->data[i + j * ptr_l->m] = ptr_a->data[i + j * ptr_a->m];
		}
	}

	// Fill P
	// Start with diagonal
	for (i=0; i<ptr_p->m; i++) ptr_p->data[i + i * ptr_p->m] = 1.0;
	// Swap the lines in P
	// indexes stored in piv_array start at 1, fixed by subtracting 1 
	double *tempvec = (double*) malloc(ptr_p->n * sizeof(double));
	for (i = min_mn - 1; i>= 0; i--) {
		// exchange  rows i and piv_array[i]
		if (i != piv_array[i]-1) {
			for (j = 0; j<ptr_p->n; j++) {
				tempvec[j] = ptr_p->data[i + j * ptr_p->m];
				ptr_p->data[i + j * ptr_p->m] = ptr_p->data[piv_array[i]-1 + j * ptr_p->m];
				ptr_p->data[piv_array[i]-1 + j * ptr_p->m] = tempvec[j];
			}
		}
	}
	free(tempvec);

	// Free
	free(piv_array);

	// Build array with results
	VALUE array = rb_ary_new2(3);
	rb_ary_store(array, 0, p);
	rb_ary_store(array, 1, l);
	rb_ary_store(array, 2, u);

	return array;
}
#endif

#ifdef LAPACK_FOUND
VALUE rb_cMatrix_inverse(VALUE self) {
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Check if is square
	matrix_raise_if_not_square(ptr);
	// Matrix a, contains input
	VALUE a = rb_cMatrix_deep_copy(self);
	MatrixData *ptr_a;
	Data_Get_Struct(a,MatrixData, ptr_a);
	int info[1];
	int* piv_array;
	piv_array = (int*) malloc(ptr_a->n * sizeof(int));
	dgetrf_(&(ptr_a->m), &(ptr_a->n), ptr_a->data, &(ptr_a->m), piv_array, info);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack LU decomposition failed");

	// Check for near-singularity (U diagonal element < threshold)
	int i;
	for (i=0; i<ptr->m; i++) {
		// U diagonal
		if (fabs(ptr_a->data[i + i * ptr_a->m]) < 1.0e-7) {
			rb_raise(rb_eRuntimeError, "Matrix is nearly singular, solution numerically unstable");
		}
	}

	int workspace_size = 1;
	dbl *workspace;
	workspace = (dbl*) malloc(workspace_size * sizeof(dbl));
	int minus_one = -1;

	// Query
	dgetri_(&(ptr_a->n), ptr_a->data, &(ptr_a->n), piv_array, workspace, &minus_one, info);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack query for workspace for matrix inverse failed");
	workspace_size = (int)*workspace; // size of new workspace
	free(workspace);
	workspace = (dbl*) malloc(workspace_size * sizeof(dbl));


	// Calculation
	dgetri_(&(ptr_a->n), ptr_a->data, &(ptr_a->n), piv_array, workspace, &workspace_size, info);
	free(piv_array);
	free(workspace);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack inverse failed");

	return a;
}
#endif

#ifdef LAPACK_FOUND
VALUE rb_cMatrix_inverse_inplace(VALUE self) {
	// Get matrix data
	MatrixData *ptr_a;
	Data_Get_Struct(self,MatrixData, ptr_a);
	// Check if is square
	matrix_raise_if_not_square(ptr_a);
	int info[1];
	int* piv_array;
	piv_array = (int*) malloc(ptr_a->n * sizeof(int));
	dgetrf_(&(ptr_a->m), &(ptr_a->n), ptr_a->data, &(ptr_a->m), piv_array, info);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack LU decomposition failed");

	int workspace_size = 1;
	dbl *workspace;
	workspace = (dbl*) malloc(workspace_size * sizeof(dbl));
	int minus_one = -1;

	// Query
	dgetri_(&(ptr_a->n), ptr_a->data, &(ptr_a->n), piv_array, workspace, &minus_one, info);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack query for workspace for matrix inverse failed");
	workspace_size = (int)*workspace; // size of new workspace
	free(workspace);
	workspace = (dbl*) malloc(workspace_size * sizeof(dbl));


	// Calculation
	dgetri_(&(ptr_a->n), ptr_a->data, &(ptr_a->n), piv_array, workspace, &workspace_size, info);
	free(piv_array);
	free(workspace);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack inverse failed");

	return Qnil;
}
#endif

#ifdef LAPACK_FOUND
VALUE rb_cMatrix_determinant(VALUE self) {
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Check if is square
	matrix_raise_if_not_square(ptr);
	// Matrix a, contains input
	VALUE a = rb_cMatrix_deep_copy(self);
	MatrixData *ptr_a;
	Data_Get_Struct(a,MatrixData, ptr_a);
	int info[1];
	int* piv_array;
	piv_array = (int*) malloc(ptr_a->n * sizeof(int));
	dgetrf_(&(ptr_a->m), &(ptr_a->n), ptr_a->data, &(ptr_a->m), piv_array, info);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack LU decomposition failed");
	double det = -1.0;
	int i;
	for (i=0; i<ptr_a->m; i++) {
		det *= ptr_a->data[i + i * ptr_a->m];
		if (piv_array[i] == i) det *= -1.0;
	}

	free(piv_array);

	return rb_float_new(det);
}
#endif

#ifdef LAPACK_FOUND
VALUE rb_cMatrix_qr2(VALUE self) {
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Matrix a, contains input
	VALUE a = rb_cMatrix_deep_copy(self);
	MatrixData *ptr_a;
	Data_Get_Struct(a,MatrixData, ptr_a);

	// Calculation
	// LAPACK directly: DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
	
	int min_mn = ptr_a->m; if (ptr_a->n < ptr_a->m) min_mn = ptr_a->n;
	double *tau_array = (double*) malloc(min_mn * sizeof(double));

	int info[1];
	int workspace_size = 1;
	dbl *workspace;
	workspace = (dbl*) malloc(workspace_size * sizeof(dbl));
	int minus_one = -1;

	dgeqrf_(&(ptr_a->m), &(ptr_a->n), ptr_a->data, &(ptr_a->m), tau_array, workspace, &minus_one, info);
	workspace_size = (int)*workspace; // size of new workspace
	free(workspace);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack query for QR workspace failed");

	workspace = (dbl*) malloc(workspace_size * sizeof(dbl));
	dgeqrf_(&(ptr_a->m), &(ptr_a->n), ptr_a->data, &(ptr_a->m), tau_array, workspace, &workspace_size, info);
	free(workspace);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack QR calculation failed");

	// Allocate output matrix Q
	VALUE q = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_q;
	Data_Get_Struct(q,MatrixData, ptr_q);
	MatrixData_alloc_zero(ptr_q, ptr_a->m, ptr_a->m);

	// Allocate output matrix R
	VALUE r = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_r;
	Data_Get_Struct(r,MatrixData, ptr_r);
	MatrixData_alloc_zero(ptr_r, ptr_a->m, ptr_a->n);

	// Copy R, prepare vectors in Q
	int i,j;
	for (i=0; i<ptr_a->m; i++) {
		for (j=0; j<ptr_a->n; j++) {
			if (j >= i) ptr_r->data[i + j * ptr_r->m] = ptr_a->data[i + j * ptr_a->m];
			if (j == i) ptr_q->data[i + j * ptr_q->m] = 1.0;
			if (j  < i) ptr_q->data[i + j * ptr_q->m] = ptr_a->data[i + j * ptr_a->m];
		}
	}
	// Multiply vectors in Q by tau^0.5
	double sqrtau;
	for (j=0; j<min_mn; j++) {
		sqrtau = sqrt(tau_array[j]);
		for (i=0; i<ptr_q->m; i++) {
			ptr_q->data[i + j * ptr_q->m] = ptr_q->data[i + j * ptr_q->m] * sqrtau;
		}
	}

	// Free
	free(tau_array);

	// Build array with results
	VALUE array = rb_ary_new2(2);
	rb_ary_store(array, 0, q);
	rb_ary_store(array, 1, r);

	return array;
}
#endif

#ifdef LAPACK_FOUND
VALUE rb_cMatrix_svd(VALUE self) {
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Matrix a, contains input
	VALUE a = rb_cMatrix_deep_copy(self);
	MatrixData *ptr_a;
	Data_Get_Struct(a,MatrixData, ptr_a);
	int m = ptr_a->m;
	int n = ptr_a->n;

	// Allocate output matrix U
	VALUE u = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_u;
	Data_Get_Struct(u,MatrixData, ptr_u);
	MatrixData_alloc(ptr_u, m, m);

	// Allocate output matrix VT
	VALUE vt = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_vt;
	Data_Get_Struct(vt ,MatrixData, ptr_vt);
	MatrixData_alloc(ptr_vt, n, n);

	int min_mn = m; if (n < m) min_mn = n;
	double *s_array = (double*) malloc(min_mn * sizeof(double));

	int info[1];
	int workspace_size = 1;
	dbl *workspace;
	workspace = (dbl*) malloc(workspace_size * sizeof(dbl));
	int minus_one = -1;

	// Query
	dgesvd_("A", "A", &m, &n, ptr_a->data, &m, s_array, ptr_u->data, &m, ptr_vt->data, &n, workspace, &minus_one, info);
	if (*info != 0) rb_raise(rb_eRuntimeError, "Lapack query for SVD workspace failed");
	workspace_size = (int)*workspace; // size of new workspace
	free(workspace);
	workspace = (dbl*) malloc(workspace_size * sizeof(dbl));

	// Calculation
	dgesvd_("A", "A", &m, &n, ptr_a->data, &m, s_array, ptr_u->data, &m, ptr_vt->data, &n, workspace, &workspace_size, info);
	free(workspace);

	// Build S matrix
	VALUE s = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_s;
	Data_Get_Struct(s,MatrixData, ptr_s);
	MatrixData_alloc_zero(ptr_s, m, n);
	int i;
	for (i = 0; i < min_mn; i++) {
		ptr_s->data[i + i * ptr_s->m] = s_array[i];
	}

	free(s_array);

	// Build final array
	VALUE array = rb_ary_new2(3);
	rb_ary_store(array, 0, u);
	rb_ary_store(array, 1, s);
	rb_ary_store(array, 2, vt);

	return array;
}
#endif

#ifdef LAPACK_FOUND
VALUE rb_cMatrix_cholesky(VALUE self) {
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Matrix a, contains input
	VALUE a = rb_cMatrix_deep_copy(self);
	MatrixData *ptr_a;
	Data_Get_Struct(a,MatrixData, ptr_a);

	// Calculation
	// LAPACK: DPOTRF( UPLO, N, A, LDA, INFO )
	int info[1];
	dpotrf_("L", &(ptr_a->n), ptr_a->data, &(ptr_a->m), info);
	if (*info != 0) rb_raise(rb_eRuntimeError, "LAPACK cholesky decomposition failed");

	// Erase the area above diagonal
	int i,j;
	for (i = 0; i < ptr_a->m; i++) {
		for (j = i+1; j < ptr_a->n; j++) {
			ptr_a->data[i + j * ptr_a->m] = 0.0;
		}
	}
	
	return a;
}
#endif

VALUE rb_cMatrix_trace(VALUE self) {
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Check if is square
	matrix_raise_if_not_square(ptr);
	// Calculate
	int i;
	double sum = 0.0;
	for (i=0; i<ptr->m; i++) sum += ptr->data[i + i*ptr->m];
	return rb_float_new(sum);
}

//==============================================================================
// Marshalling
//==============================================================================

VALUE rb_cMatrix_marshal_dump(VALUE self, VALUE depth) {
	// Get matrix data
	MatrixData *ptr;
	Data_Get_Struct(self,MatrixData, ptr);
	// Create string
	VALUE str = rb_str_buf_new(2*sizeof(int) + ptr->m * ptr->n * sizeof(dbl));
	// Wrire string
	rb_str_buf_cat(str, (const char*)(&(ptr->m)), sizeof(int)) ;
	rb_str_buf_cat(str, (const char*)(&(ptr->n)), sizeof(int)) ;
	rb_str_buf_cat(str, (const char*)(ptr->data), ptr->m * ptr->n * sizeof(dbl));
	return str;
}

VALUE rb_cMatrix_marshal_load(VALUE class, VALUE str) {
	// Convert value to C string
	char* p = (char*)RSTRING_PTR(str);
	// Get matrix size
	int* mn = (int*)p;
	// Allocate result
	VALUE result = rb_cMatrix_new(class);
	MatrixData *ptr_r;
	Data_Get_Struct(result,MatrixData, ptr_r);
	MatrixData_alloc(ptr_r, mn[0], mn[1]);
	// Copy data
	memcpy(ptr_r->data, p + 2*sizeof(int), ptr_r->m *  ptr_r->n * sizeof(dbl));
	return result;
}

//##############################################################################
// Extension init
//##############################################################################

void Init_matrix_c () {

	rb_cMatrix = rb_define_class_under(rb_mAlgebra, "Matrix", rb_cObject);

	rb_define_alloc_func(rb_cMatrix, rb_cMatrix_new);

	// Constructors
	rb_define_method(rb_cMatrix, "initialize", rb_cMatrix_initialize, 2);
	rb_define_singleton_method(rb_cMatrix, "from_array", rb_cMatrix_from_array, 1);
	rb_define_singleton_method(rb_cMatrix, "filled", rb_cMatrix_filled, 3);
	rb_define_singleton_method(rb_cMatrix, "zero", rb_cMatrix_zero, -1);
	/// zero(m,n)
	rb_define_singleton_method(rb_cMatrix, "diagonal1", rb_cMatrix_diagonal_1, 2);
	rb_define_singleton_method(rb_cMatrix, "diagonal2", rb_cMatrix_diagonal_2, 1);

	// Matrix information
	rb_define_method(rb_cMatrix,"m", rb_cMatrix_size_m, 0);
	rb_define_method(rb_cMatrix,"n", rb_cMatrix_size_n, 0);

	// Elemet access
	rb_define_method(rb_cMatrix,"[]", rb_cMatrix_get_ij, 2);
	rb_define_method(rb_cMatrix,"[]=", rb_cMatrix_set_ij, 3);
	rb_define_method(rb_cMatrix,"submatrix", rb_cMatrix_submatrix, 4);

	// Row & col access
	rb_define_method(rb_cMatrix, "row_as_array", rb_cMatrix_row_to_a, 1);
	rb_define_method(rb_cMatrix, "col_as_array", rb_cMatrix_col_to_a, 1);
	rb_define_method(rb_cMatrix, "column_as_array", rb_cMatrix_col_to_a, 1);
	rb_define_method(rb_cMatrix, "row_as_vector", rb_cMatrix_row_to_vector, 1);
	rb_define_method(rb_cMatrix, "col_as_vector", rb_cMatrix_col_to_vector, 1);
	rb_define_method(rb_cMatrix, "column_as_vector", rb_cMatrix_col_to_vector, 1);
	rb_define_method(rb_cMatrix, "diagonal_as_vector", rb_cMatrix_diagonal_to_vector, 0);

	// Misc
	//rb_define_method(rb_cMatrix, "to_s", rb_cMatrix_to_s, 0);
	rb_define_method(rb_cMatrix, "to_a", rb_cMatrix_to_a, 0);
	rb_define_method(rb_cMatrix, "deep_copy", rb_cMatrix_deep_copy, 0);
	rb_define_method(rb_cMatrix, "to_vector", rb_cMatrix_to_vector, 0);

	// Operators
	rb_define_method(rb_cMatrix, "matrix_multiply", rb_cMatrix_matrix_multiply, 1);
	rb_define_method(rb_cMatrix, "multiply_yield_diagonal", rb_cMatrix_matrix_multiply_yield_diagonal, 1);
	rb_define_method(rb_cMatrix, "matrix_vector_multiply", rb_cMatrix_matrix_vector_multiply, 1);
	rb_define_method(rb_cMatrix, "+", rb_cMatrix_plus, 1);
	rb_define_method(rb_cMatrix, "-", rb_cMatrix_minus, 1);
	rb_define_method(rb_cMatrix, "plus!", rb_cMatrix_plus_inplace, 1);
	rb_define_method(rb_cMatrix, "minus!", rb_cMatrix_minus_inplace, 1);
	rb_define_method(rb_cMatrix, "comutative_multiply", rb_cMatrix_comutative_multiply, 1);
	rb_define_method(rb_cMatrix, "mul!", rb_cMatrix_multiply_inplace, 1);

	// Comparison
	rb_define_method(rb_cMatrix, "==", rb_cMatrix_equal, 1);
	rb_define_method(rb_cMatrix, "=~", rb_cMatrix_equal_eps, 1);

	// Matrix operations
	rb_define_method(rb_cMatrix, "transpose", rb_cMatrix_transpose, 0);
	rb_define_method(rb_cMatrix, "transpose!", rb_cMatrix_transpose_inplace, 0);
	rb_define_method(rb_cMatrix, "trace", rb_cMatrix_trace, 0);

	// Matrix operations calling LAPACK
#ifdef LAPACK_FOUND
	rb_define_method(rb_cMatrix, "eigensystem", rb_cMatrix_eigensystem, 0);
	rb_define_method(rb_cMatrix, "lu", rb_cMatrix_lu, 0);
	rb_define_method(rb_cMatrix, "inverse", rb_cMatrix_inverse, 0);
	rb_define_method(rb_cMatrix, "inverse!", rb_cMatrix_inverse_inplace, 0);
	rb_define_method(rb_cMatrix, "determinant", rb_cMatrix_determinant, 0);
	rb_define_method(rb_cMatrix, "qr2", rb_cMatrix_qr2, 0);
	rb_define_method(rb_cMatrix, "svd", rb_cMatrix_svd, 0);
	rb_define_method(rb_cMatrix, "cholesky", rb_cMatrix_cholesky, 0);
#endif

	// Marshaling
	rb_define_method(rb_cMatrix, "_dump", RUBY_METHOD_FUNC(rb_cMatrix_marshal_dump), 1);
	rb_define_singleton_method(rb_cMatrix, "_load", RUBY_METHOD_FUNC(rb_cMatrix_marshal_load), 1);
}
