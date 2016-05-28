#include <ruby.h>
#include <math.h>
#include <stdio.h>

#ifdef UMFPACK_FOUND
#include <umfpack.h>
#endif

#include "algebra_c.h"
#include "sparse_matrix_csx_c.h"

//##############################################################################
// Sparse matrix data structure
//##############################################################################

SparseMatrixCSXData *SparseMatrixCSXData_new() {
	SparseMatrixCSXData *ptr = (SparseMatrixCSXData*) malloc(sizeof(SparseMatrixCSXData));
	ptr->row_mode = 0; // 1 for row, 0 for col
	ptr->m = 0;
	ptr->n = 0;
	ptr->nnz = 0;
	ptr->data_values = NULL;
	ptr->data_records = NULL;
	ptr->data_indices = NULL;
	return ptr;
}

void SparseMatrixCSXData_free(SparseMatrixCSXData* md) {
	if (md->data_values != NULL) free(md->data_values);
	if (md->data_records != NULL) free(md->data_records);
	if (md->data_indices != NULL) free(md->data_indices);
	free(md);
}

void SparseMatrixCSXData_alloc(SparseMatrixCSXData *ptr, int row_mode, int m, int n, int nnz, int recs) {
	ptr->row_mode = row_mode;
	ptr->m = m;
	ptr->n = n;
	ptr->nnz = nnz;
	ptr->data_values = (dbl*) malloc(nnz * sizeof(dbl));
	ptr->data_indices = (int*) malloc(nnz * sizeof(int));
	ptr->data_records = (int*) malloc((recs+1) * sizeof(int));
}

//##############################################################################
// Helpers
//##############################################################################

void sparse_matrix_raise_index_error(SparseMatrixCSXData *a, int i, int j) {
	if (i >= a->m || i < 0) {
		rb_raise(rb_eIndexError, "Row index (i) out of bounds");
	}
	if (j >= a->n || j < 0) {
		rb_raise(rb_eIndexError, "Column index (j) out of bounds");
	}
}

// Classes defined here
VALUE rb_cSparseMatrixCSX;

//##############################################################################
// SparseMatrix class methods
//##############################################################################

//==============================================================================
// Allocation function
//==============================================================================

VALUE rb_cSparseMatrixCSX_new(VALUE class) {
	SparseMatrixCSXData *ptr = SparseMatrixCSXData_new();
	VALUE data = Data_Wrap_Struct(class, 0, SparseMatrixCSXData_free, ptr);
	return data;
}

//==============================================================================
// Constructors
//==============================================================================

VALUE rb_cSparseMatrixCSX_initialize(VALUE self, VALUE mode, VALUE colsize, VALUE rowsize, VALUE vals, VALUE indices, VALUE record_i) {
	SparseMatrixCSXData *ptr;
	int m = NUM2INT(colsize);
	int n = NUM2INT(rowsize);
	int nnz = RARRAY_LEN(vals);
	int row_mode;
	int i;
	int recs;

	// Set mode
	if (mode == ID2SYM(rb_intern("row"))) {
		row_mode = 1;
	} else {
		if ID2SYM((mode == rb_intern("col"))) {
			row_mode = 0;
		} else {
			rb_raise(rb_eSparseTypeError, "Row mode must be either :row or :col");
		}
	}	

	// Check dimensions of the arrays
	if (row_mode) {
		if(RARRAY_LEN(record_i) != m+1) rb_raise(rb_eRuntimeError, "Size of row indices array must be m+1");
		if(RARRAY_LEN(vals) != RARRAY_LEN(indices)) rb_raise(rb_eRuntimeError, "Size of values and column indices arrays must be the same");
		// Number of records = number of rows
		recs = m;
	} else {
		if(RARRAY_LEN(record_i) != n+1) rb_raise(rb_eRuntimeError, "Size of col indices array must be n+1");
		if(RARRAY_LEN(vals) != RARRAY_LEN(indices)) rb_raise(rb_eRuntimeError, "Size of values and row indices arrays must be the same");
		// Number of records = number of cols
		recs = n;
	}


	// Allocate arrays
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	SparseMatrixCSXData_alloc(ptr, row_mode, m, n, nnz, recs);

	// Copy array of values
	VALUE* array_ptr = RARRAY_PTR(vals);
	for(i = 0; i < nnz; i++) {
		ptr->data_values[i] = NUM2DBL(array_ptr[i]);
	}

	// Copy array of data indexes
	VALUE* array_ptr_c = RARRAY_PTR(indices);
	for(i = 0; i < nnz; i++) {
		ptr->data_indices[i] = NUM2INT(array_ptr_c[i]);
	}

	// Copy array of record indexes
	VALUE* array_ptr_r = RARRAY_PTR(record_i);
	for(i = 0; i < (recs+1); i++) {
		ptr->data_records[i] = NUM2INT(array_ptr_r[i]);
	}
	return Qnil;
}

//==============================================================================
// Basics
//==============================================================================

VALUE rb_cSparseMatrixCSX_deep_copy(VALUE self) {
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	int recs;
	if (ptr->row_mode) {
		// Number of records = number of rows
		recs = ptr->m;
	} else {
		// Number of records = number of cols
		recs = ptr->n;
	}
	// Allocate result
	VALUE result = rb_cSparseMatrixCSX_new(rb_class_of(self));
	SparseMatrixCSXData *ptr_r;
	Data_Get_Struct(result, SparseMatrixCSXData, ptr_r);
	SparseMatrixCSXData_alloc(ptr_r, ptr->row_mode, ptr->m, ptr->n, ptr->nnz, recs);
	// Memcpy on the arrays
	memcpy(ptr_r->data_values, ptr->data_values, ptr->nnz * sizeof(dbl));
	memcpy(ptr_r->data_indices, ptr->data_indices, ptr->nnz * sizeof(int));
	memcpy(ptr_r->data_records, ptr->data_records, (recs + 1) * sizeof(int));
	return result;
}

VALUE rb_cSparseMatrixCSX_rewrite(VALUE self, VALUE mode) {
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	int newmode, newrecs;
	if (mode == ID2SYM(rb_intern("row"))) {
		newmode = 1;
		newrecs = ptr->m;
	} else {
		if ID2SYM((mode == rb_intern("col"))) {
			newmode = 0;
			newrecs = ptr->n;
		} else {
			rb_raise(rb_eSparseTypeError, "Row mode must be either :row or :col");
		}
	}
	// No mode change - use deep copy
	if (newmode == ptr->row_mode) {
		return rb_cSparseMatrixCSX_deep_copy(self);
	}

	// Allocate result
	VALUE result = rb_cSparseMatrixCSX_new(rb_class_of(self));
	SparseMatrixCSXData *ptr_r;
	Data_Get_Struct(result, SparseMatrixCSXData, ptr_r);
	SparseMatrixCSXData_alloc(ptr_r, newmode, ptr->m, ptr->n, ptr->nnz, newrecs);

	// Traverse nonzero elements and write new data
	int count;
	int reci;
	int i, j, c, rs, re;
	int size1, size2;
	if (newmode) {
		size1 = ptr->m;
		size2 = ptr->n;
	} else {
		size1 = ptr->n;
		size2 = ptr->m;
	}
	count = 0;
	reci = 0;
	rb_funcall(self, rb_intern("raise_slow"), 0);
	for(i=0; i<size1; i++) {
		ptr_r->data_records[reci] = count;
		reci++;
		// Each nonzero with index:
		for (j = 0; j < size2; j++) {
			// Record start and end
			rs = ptr->data_records[j];
			re = ptr->data_records[j+1];
			// Traverse the column
			for(c = rs; c < re; c++) {
				if(i ==  ptr->data_indices[c]) {
					ptr_r->data_values[count] = ptr->data_values[c];
					ptr_r->data_indices[count] = j;
					count++;
					break;
				}
			}
		}
	}
	ptr_r->data_records[reci] = count;
	
	// Result
	return result;
}

//==============================================================================
// Matrix information
//==============================================================================

VALUE rb_cSparseMatrixCSX_type(VALUE self) {
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	if (ptr->row_mode) {
		return ID2SYM(rb_intern("row"));
	} else {
		return ID2SYM(rb_intern("col"));
	}
}

VALUE rb_cSparseMatrixCSX_change_type_inplace(VALUE self, VALUE mode) {
	int oldmode, tmp;
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	oldmode = ptr->row_mode;
	// Set mode
	if (mode == ID2SYM(rb_intern("row"))) {
		ptr->row_mode = 1;
	} else {
		if ID2SYM((mode == rb_intern("col"))) {
			ptr->row_mode = 0;
		} else {
			rb_raise(rb_eSparseTypeError, "Row mode must be either :row or :col");
		}
	}
	// If mode changes
	if (ptr->row_mode != oldmode) {
		// swap m and n
		tmp = ptr->m;
		ptr->m = ptr->n;
		ptr->n = tmp;
	}
	return Qnil;
}

VALUE rb_cSparseMatrixCSX_size_m(VALUE self) {
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	return INT2NUM(ptr->m);
}

VALUE rb_cSparseMatrixCSX_size_n(VALUE self) {
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	return INT2NUM(ptr->n);
}

//==============================================================================
// Element access
//==============================================================================

VALUE rb_cSparseMatrixCSX_get_ij(VALUE self, VALUE i, VALUE j) {
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	int ii = NUM2INT(i);
	int jj = NUM2INT(j);
	int ri, rj;
	// Out of bounds check
	sparse_matrix_raise_index_error(ptr, ii, jj);

	// row/col mode
	if (ptr->row_mode) {
		ri = ii;
		rj = jj;
	} else {
		ri = jj;
		rj = ii;
	}

	int rs, re, c;
	// Record start and end
	rs = ptr->data_records[ri];
	re = ptr->data_records[ri+1];
	// Traverse the row
	for(c = rs; c < re; c++) {
		if(rj ==  ptr->data_indices[c]) return rb_float_new(ptr->data_values[c]);
	}

	// Element not found -> return zero
	return rb_float_new(0.0);
}

VALUE rb_cSparseMatrixCSX_set_ij(VALUE self, VALUE i, VALUE j, VALUE v) {
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	int ii = NUM2INT(i);
	int jj = NUM2INT(j);
	int ri, rj;
	// Out of bounds check
	sparse_matrix_raise_index_error(ptr, ii, jj);

	// row/col mode
	if (ptr->row_mode) {
		ri = ii;
		rj = jj;
	} else {
		ri = jj;
		rj = ii;
	}

	int rs, re, c;
	// Record start and end
	rs = ptr->data_records[ri];
	re = ptr->data_records[ri+1];
	// Traverse the row
	for(c = rs; c < re; c++) {
		if(rj == ptr->data_indices[c]) {
			ptr->data_values[c] = NUM2DBL(v);
			return v;
		}
	}

	// Element not found -> raise error
	rb_raise(rb_eSparseWriteError, "SparseMatrix can not write into element which is zero");
}

//==============================================================================
// Row / col access
//==============================================================================

VALUE rb_cSparseMatrixCSX_row_to_a(VALUE self, VALUE index) {
	int i, j, c, rs, re;
	// Get matrix data
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	// Check bounds
	i = NUM2INT(index);
	sparse_matrix_raise_index_error(ptr, i, 0);

	// New array
	VALUE array = rb_ary_new2(ptr->n);
	// Fill array with zeroes
	for (j = 0; j < ptr->n; j++) {
		rb_ary_store(array, j, rb_float_new(0.0)); 
	}

	// row/col mode
	if (ptr->row_mode) {
		// Row start and end
		rs = ptr->data_records[i];
		re = ptr->data_records[i+1];
		// Traverse the row and write values
		for(j = rs; j < re; j++) {
			rb_ary_store(array, ptr->data_indices[j], rb_float_new(ptr->data_values[j]));
		}
	} else {
		// Get and write values
		rb_funcall(self, rb_intern("raise_slow"), 0);
		for (j = 0; j < ptr->n; j++) {
			// Record (col) start and end
			rs = ptr->data_records[j];
			re = ptr->data_records[j+1];
			// Traverse the column
			for(c = rs; c < re; c++) {
				if(i ==  ptr->data_indices[c]) {
					 rb_ary_store(array, j, rb_float_new(ptr->data_values[c]));
					 break;
				}
			}
		}
	}
	return array;
}

VALUE rb_cSparseMatrixCSX_col_to_a(VALUE self, VALUE index) {
	int i, j, c, rs, re;
	// Get matrix data
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	// Check bounds
	i = NUM2INT(index);
	sparse_matrix_raise_index_error(ptr, 0, i);

	// New array
	VALUE array = rb_ary_new2(ptr->m);
	// Fill array with zeroes
	for (j = 0; j < ptr->m; j++) {
		rb_ary_store(array, j, rb_float_new(0.0)); 
	}

	// row/col mode
	if (ptr->row_mode) {
		// Slow search for each element
		rb_funcall(self, rb_intern("raise_slow"), 0);
		for (j = 0; j < ptr->m; j++) {
			// Record (col) start and end
			rs = ptr->data_records[j];
			re = ptr->data_records[j+1];
			// Traverse the column
			for(c = rs; c < re; c++) {
				if(i ==  ptr->data_indices[c]) {
					rb_ary_store(array, j, rb_float_new(ptr->data_values[c]));
					break;
				}
			}
		}
	} else {
		// Row start and end
		rs = ptr->data_records[i];
		re = ptr->data_records[i+1];
		// Traverse the row and write values
		for(j = rs; j < re; j++) {
			rb_ary_store(array, ptr->data_indices[j], rb_float_new(ptr->data_values[j]));
		}
	}
	return array;
}

//==============================================================================
// Iterators
//==============================================================================

VALUE rb_cSparseMatrixCSX_row_each_nonzero_with_index(VALUE self, VALUE index) {
	int i, j, c, rs, re;
	// Get matrix data
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	// Check bounds
	i = NUM2INT(index);
	sparse_matrix_raise_index_error(ptr, i, 0);

	// New array which is yielded - size 2
	VALUE array = rb_ary_new2(2);

	// row/col mode
	if (ptr->row_mode) {
		// Row start and end
		rs = ptr->data_records[i];
		re = ptr->data_records[i+1];
		// Traverse the row and write values
		for(j = rs; j < re; j++) {
			rb_ary_store(array, 0, rb_float_new(ptr->data_values[j])); 
			rb_ary_store(array, 1, INT2NUM(ptr->data_indices[j]));
			rb_yield(array);
		}
	} else {
		// Slow search for each element
		rb_funcall(self, rb_intern("raise_slow"), 0);
		for (j = 0; j < ptr->n; j++) {
			// Record (col) start and end
			rs = ptr->data_records[j];
			re = ptr->data_records[j+1];
			// Traverse the column
			for(c = rs; c < re; c++) {
				if(i ==  ptr->data_indices[c]) {
					rb_ary_store(array, 0, rb_float_new(ptr->data_values[c])); 
					rb_ary_store(array, 1, INT2NUM(j));
					rb_yield(array);
					break;
				}
			}
		}
	}
	return Qnil;
}

VALUE rb_cSparseMatrixCSX_col_each_nonzero_with_index(VALUE self, VALUE index) {
	int i, j, c, rs, re;
	// Get matrix data
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	// Check bounds
	i = NUM2INT(index);
	sparse_matrix_raise_index_error(ptr, 0, i);

	// New array which is yielded - size 2
	VALUE array = rb_ary_new2(2);

	// row/col mode
	if (ptr->row_mode) {
		// Slow search for each element
		rb_funcall(self, rb_intern("raise_slow"), 0);
		for (j = 0; j < ptr->m; j++) {
			// Record (col) start and end
			rs = ptr->data_records[j];
			re = ptr->data_records[j+1];
			// Traverse the column
			for(c = rs; c < re; c++) {
				if(i ==  ptr->data_indices[c]) {
					rb_ary_store(array, 0, rb_float_new(ptr->data_values[c])); 
					rb_ary_store(array, 1, INT2NUM(j));
					rb_yield(array);
					break;
				}
			}
		}
	} else {
		// Row start and end
		rs = ptr->data_records[i];
		re = ptr->data_records[i+1];
		// Traverse the row and write values
		for(j = rs; j < re; j++) {
			rb_ary_store(array, 0, rb_float_new(ptr->data_values[j])); 
			rb_ary_store(array, 1, INT2NUM(ptr->data_indices[j]));
			rb_yield(array);
		}
	}
	return Qnil;
}

//==============================================================================
// In-matrix operations
//==============================================================================

VALUE rb_cSparseMatrixCSX_record_dot_record(VALUE self, VALUE index1, VALUE index2) {
	// Get matrix data
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	// Convert arguments (bounds check done in ruby, this method is not called directly)
	int i1, i2;
	i1 = NUM2INT(index1);
	i2 = NUM2INT(index2);

	// Traverse the rows simultaneously, using two pointing indices
	int re1, re2, p1, p2;
	int di1, di2;

	p1 = ptr->data_records[i1];
	re1 = ptr->data_records[i1+1];

	p2 = ptr->data_records[i2];
	re2 = ptr->data_records[i2+1];

	dbl sum = 0.0;
	while (p1 < re1 && p2 < re2) {
		di1 = ptr->data_indices[p1];
		di2 = ptr->data_indices[p2];
		if (di1 == di2) {
			sum += ptr->data_values[p1] * ptr->data_values[p2];
			p1++;
			p2++;
		} else {
			if (di1 < di2) p1++; else p2++;
		}
	}

	return rb_float_new(sum);
}

//==============================================================================
// Access to  array representation
//==============================================================================

VALUE rb_cSparseMatrixCSX_csx_array_values(VALUE self){
	SparseMatrixCSXData *ptr;
	int i;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);

	VALUE array = rb_ary_new2(ptr->nnz);
	for (i = 0; i < ptr->nnz; i++) {
		rb_ary_store(array, i, rb_float_new(ptr->data_values[i]));
	}

	return array;
}

VALUE rb_cSparseMatrixCSX_csx_array_indices(VALUE self){
	SparseMatrixCSXData *ptr;
	int i;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);

	VALUE array = rb_ary_new2(ptr->nnz);
	for (i = 0; i < ptr->nnz; i++) {
		rb_ary_store(array, i, INT2NUM(ptr->data_indices[i]));
	}

	return array;
}

VALUE rb_cSparseMatrixCSX_csx_array_records(VALUE self){
	SparseMatrixCSXData *ptr;
	int i;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);

	VALUE array = rb_ary_new2(ptr->m+1);
	for (i = 0; i < (ptr->m+1); i++) {
		rb_ary_store(array, i, INT2NUM(ptr->data_records[i]));
	}

	return array;
}

//==============================================================================
// Type conversion
//==============================================================================

VALUE rb_cSparseMatrixCSX_to_matrix(VALUE self){
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	// Allocate result matrix
	VALUE result = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_m;
	Data_Get_Struct(result,MatrixData, ptr_m);
	MatrixData_alloc_zero(ptr_m, ptr->m, ptr->n);
	if (ptr->row_mode) {
		// Iterate
		int i,ji, rs, re;
		for (i=0; i<ptr->m; i++) {
			// Row start and end
			rs = ptr->data_records[i];
			re = ptr->data_records[i+1];
			for(ji = rs; ji < re; ji++) {
				// index ptr->data_indices[ji]
				// value ptr->data_values[ji]
				ptr_m->data[i + ptr->data_indices[ji] * ptr_m->m] = ptr->data_values[ji];
			}
		}
		return result;
	} else {
		// Iterate
		int j, ii, cs, ce;
		for (j=0; j<ptr->n; j++) {
			// Col start and end
			cs = ptr->data_records[j];
			ce = ptr->data_records[j+1];
			for(ii = cs; ii < ce; ii++) {
				// index ptr->data_indices[ii]
				// value ptr->data_values[ii]
				ptr_m->data[ptr->data_indices[ii] + j * ptr_m->m] = ptr->data_values[ii];
			}
		}
		return result;
	}
}

VALUE rb_cSparseMatrixCSX_to_sparse_matrix(VALUE self){
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	// Allocate result matrix
	VALUE result;
	result = rb_funcall(rb_eval_string("SparseMatrix"), rb_intern("new"), 2, INT2NUM(ptr->m), INT2NUM(ptr->n));
	// Assignment method
	VALUE assignment = rb_intern("[]=");
	if (ptr->row_mode) {
		// Iterate - easy
		int i,ji, rs, re;
		for (i=0; i<ptr->m; i++) {
			// Row start and end
			rs = ptr->data_records[i];
			re = ptr->data_records[i+1];
			for(ji = rs; ji < re; ji++) {
				// index ptr->data_indices[ji]
				// value ptr->data_values[ji]
				rb_funcall(result, assignment, 3, INT2NUM(i),  INT2NUM(ptr->data_indices[ji]), rb_float_new(ptr->data_values[ji]));
			}
		}
		return result;
	} else {
		// Iterate
		int j, ii, cs, ce;
		for (j=0; j<ptr->n; j++) {
			// Col start and end
			cs = ptr->data_records[j];
			ce = ptr->data_records[j+1];
			for(ii = cs; ii < ce; ii++) {
				// index ptr->data_indices[ii]
				// value ptr->data_values[ii]
				rb_funcall(result, assignment, 3, INT2NUM(ptr->data_indices[ii]), INT2NUM(j), rb_float_new(ptr->data_values[ii]));
			}
		}
		return result;
	}
}

//==============================================================================
// Operators
//==============================================================================

VALUE rb_cSparseMatrixCSX_matrix_vector_multiply_c(VALUE self, VALUE vector) {
	vector_raise_if_not_vector(vector);
	// Get matrix data
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	if (ptr->row_mode) {
		// Get vector data
		VectorData *ptr_v;
		Data_Get_Struct(vector,VectorData, ptr_v);
		// Check size
		if (ptr->n != ptr_v->size) rb_raise(rb_eDimensionError, "Cant multiply matrix with vector of wrong size");
		// Result variables
		VALUE result;
		VectorData *ptr_r;
		int i,ji, rs, re;
		double sum;
		// The result is Vector
		result = rb_cVector_new(rb_class_of(vector));
		Data_Get_Struct(result,VectorData, ptr_r);
		VectorData_alloc(ptr_r, ptr->m);
		// Multiply
		for (i=0; i<ptr->m; i++) {
			sum = 0.0;
			// Row start and end
			rs = ptr->data_records[i];
			re = ptr->data_records[i+1];
			for(ji = rs; ji < re; ji++) {
				// index ptr->data_indices[ji]
				// value ptr->data_values[ji]
				sum += ptr_v->data[ptr->data_indices[ji]] * ptr->data_values[ji];
			}
			ptr_r->data[i] = sum;
		}
		return result;
	} else {
		rb_raise(rb_eRuntimeError, "CSC matrix multiplication not implemented");
	}
}

//==============================================================================
// UMFPACK solver
//==============================================================================

#ifdef UMFPACK_FOUND
VALUE rb_cSparseMatrixCSX_solve_umfpack(VALUE self, VALUE vector){
	SparseMatrixCSXData *ptr;
	Data_Get_Struct(self,SparseMatrixCSXData, ptr);
	VectorData *ptr_v;
	Data_Get_Struct(vector,VectorData, ptr_v);

	// Allocate result
	VALUE result = rb_cVector_new(rb_cVector);
	VectorData *ptr_x;
	Data_Get_Struct(result,VectorData, ptr_x);
	VectorData_alloc(ptr_x, ptr_v->size);
	// Call UMFPACK routines
	if (ptr->row_mode) {
		dbl *null = (dbl*) NULL;
		void *Symbolic, *Numeric;
		int retval;
		retval = umfpack_di_symbolic (ptr->m, ptr->m, ptr->data_records, ptr->data_indices, ptr->data_values, &Symbolic, null, null);   // symbolic reordering of the matrix given by CAp, CAi, and CAx
		switch(retval){
			case UMFPACK_ERROR_n_nonpositive: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_n_nonpositive"); break;
			case UMFPACK_ERROR_invalid_matrix: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_invalid_matrix"); break;
			case UMFPACK_ERROR_out_of_memory: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_out_of_memory"); break;
			case UMFPACK_ERROR_argument_missing: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_argument_missing"); break;
			case UMFPACK_ERROR_internal_error: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_internal_error"); break;
		}
		retval = umfpack_di_numeric (ptr->data_records, ptr->data_indices, ptr->data_values, Symbolic, &Numeric, null, null);   // factorization of CAp, CAi, CAx
		switch(retval){
			case UMFPACK_WARNING_singular_matrix: rb_raise(rb_eRuntimeError, "UMFPACK_WARNING_singular_matrix"); break;
			case UMFPACK_ERROR_out_of_memory: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_out_of_memory"); break;
			case UMFPACK_ERROR_argument_missing: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_argument_missing"); break;
			case UMFPACK_ERROR_invalid_Symbolic_object: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_invalid_Symbolic_object"); break;
			case UMFPACK_ERROR_different_pattern: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_different_pattern"); break;
		}
		umfpack_di_free_symbolic (&Symbolic) ;
		retval = umfpack_di_solve (UMFPACK_A, ptr->data_records, ptr->data_indices, ptr->data_values, ptr_x->data, ptr_v->data, Numeric, null, null);
		switch(retval){
			case UMFPACK_WARNING_singular_matrix: rb_raise(rb_eRuntimeError, "UMFPACK_WARNING_singular_matrix"); break;
			case UMFPACK_ERROR_out_of_memory: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_out_of_memory"); break;
			case UMFPACK_ERROR_argument_missing: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_argument_missing"); break;
			case UMFPACK_ERROR_invalid_system: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_invalid_system"); break;
			case UMFPACK_ERROR_invalid_Numeric_object: rb_raise(rb_eRuntimeError, "UMFPACK_ERROR_invalid_Numeric_object"); break;
		}
		umfpack_di_free_numeric (&Numeric);
		return result;
	} else {
		rb_raise(rb_eRuntimeError, "CSC matrix solver not implemented");
	}

}
#endif

//##############################################################################
// Extension init
//##############################################################################

void Init_sparse_matrix_csx_c () {

	rb_cSparseMatrixCSX = rb_define_class_under(rb_mAlgebra, "SparseMatrixCSX", rb_cObject);

	rb_define_alloc_func(rb_cSparseMatrixCSX, rb_cSparseMatrixCSX_new);

	// Constructors
	rb_define_method(rb_cSparseMatrixCSX, "initialize", rb_cSparseMatrixCSX_initialize, 6);

	// Basics
	rb_define_method(rb_cSparseMatrixCSX,"deep_copy", rb_cSparseMatrixCSX_deep_copy, 0);
	rb_define_method(rb_cSparseMatrixCSX,"rewrite", rb_cSparseMatrixCSX_rewrite, 1);

	// Matrix information
	rb_define_method(rb_cSparseMatrixCSX,"type", rb_cSparseMatrixCSX_type, 0);
	rb_define_method(rb_cSparseMatrixCSX,"change_type!", rb_cSparseMatrixCSX_change_type_inplace, 1);
	rb_define_method(rb_cSparseMatrixCSX,"m", rb_cSparseMatrixCSX_size_m, 0);
	rb_define_method(rb_cSparseMatrixCSX,"n", rb_cSparseMatrixCSX_size_n, 0);

	// Elemet access
	rb_define_method(rb_cSparseMatrixCSX,"[]", rb_cSparseMatrixCSX_get_ij, 2);
	rb_define_method(rb_cSparseMatrixCSX,"[]=", rb_cSparseMatrixCSX_set_ij, 3);

	// Row & col access
	rb_define_method(rb_cSparseMatrixCSX, "row_as_array", rb_cSparseMatrixCSX_row_to_a, 1);
	rb_define_method(rb_cSparseMatrixCSX, "col_as_array", rb_cSparseMatrixCSX_col_to_a, 1);

	// Iterators
	rb_define_method(rb_cSparseMatrixCSX, "row_each_nonzero_with_index", rb_cSparseMatrixCSX_row_each_nonzero_with_index, 1);
	rb_define_method(rb_cSparseMatrixCSX, "col_each_nonzero_with_index", rb_cSparseMatrixCSX_col_each_nonzero_with_index, 1);

	// In-matrix operations
	rb_define_method(rb_cSparseMatrixCSX, "record_dot_record", rb_cSparseMatrixCSX_record_dot_record, 2);

	// Access to  arrays
	rb_define_method(rb_cSparseMatrixCSX,"csx_array_values", rb_cSparseMatrixCSX_csx_array_values, 0);
	rb_define_method(rb_cSparseMatrixCSX,"csx_array_indices", rb_cSparseMatrixCSX_csx_array_indices, 0);
	rb_define_method(rb_cSparseMatrixCSX,"csx_array_records", rb_cSparseMatrixCSX_csx_array_records, 0);

	// Type conversion
	rb_define_method(rb_cSparseMatrixCSX,"to_matrix", rb_cSparseMatrixCSX_to_matrix, 0);
	rb_define_method(rb_cSparseMatrixCSX,"to_sparse_matrix", rb_cSparseMatrixCSX_to_sparse_matrix, 0);

	// Operators
	rb_define_method(rb_cSparseMatrixCSX, "matrix_vector_multiply_c", rb_cSparseMatrixCSX_matrix_vector_multiply_c, 1);

	// UMFPACK solver
#ifdef UMFPACK_FOUND
	rb_define_method(rb_cSparseMatrixCSX,"solve", rb_cSparseMatrixCSX_solve_umfpack, 1);
#endif

}
