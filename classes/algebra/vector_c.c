#include <ruby.h>
#include <math.h>

#include "algebra_c.h"

// Classes defined here
VALUE rb_cVector;

VectorData *VectorData_new() {
	VectorData *ptr = (VectorData*) malloc(sizeof(VectorData));
	ptr->size = 0;
	ptr->data = NULL;
	return ptr;
}

void VectorData_free(VectorData* md) {
	if (md->data != NULL) free(md->data);
	free(md);
}

void VectorData_alloc(VectorData *ptr, int n) {
	ptr->size = n;
	ptr->data = (dbl*) malloc(n * sizeof(dbl));
}

void VectorData_alloc_zero(VectorData *ptr, int n) {
	VectorData_alloc(ptr,  n);
	int i;
	for(i = 0; i < n; i++) ptr->data[i] = 0.0;
}

//##############################################################################
// Helpers
//##############################################################################

void vector_raise_if_not_same_size(VectorData *a, VectorData *b) {
	if (a->size != b->size) {
		rb_raise(rb_eDimensionError, "Vectors must be of same size");
	}
}

void vector_raise_index_error(VectorData *a, int index) {
	if (index >= a->size || index < 0) {
		rb_raise(rb_eIndexError, "Index out of bounds");
	}
}

void vector_raise_if_not_vector(VALUE val) {
	if (rb_obj_is_kind_of(val, rb_cVector) != Qtrue) {
		rb_raise(rb_eTypeError, "Argument is not a Vector");
	}
}

//##############################################################################
// Vector class methods
//##############################################################################

//==============================================================================
// Allocation function
//==============================================================================

VALUE rb_cVector_new(VALUE class) {
	VectorData *ptr = VectorData_new();
	VALUE data = Data_Wrap_Struct(class, 0, VectorData_free, ptr);
	return data;
}

//==============================================================================
// Constructors
//==============================================================================

VALUE rb_cVector_initialize(VALUE self, VALUE size) {
	VectorData *ptr;
	int n = NUM2INT(size);
	if (n <= 0) rb_raise(rb_eTypeError, "Vector can not have size 0 or smaller");
	Data_Get_Struct(self,VectorData, ptr);
	VectorData_alloc_zero(ptr, n);
	return Qnil;
}

VALUE rb_cVector_zero(VALUE class, VALUE size) {
	VALUE result = rb_cVector_new(class);
	VectorData *ptr;
	int n = NUM2INT(size);
	if (n <= 0) rb_raise(rb_eTypeError, "Vector can not have size 0 or smaller");
	Data_Get_Struct(result,VectorData, ptr);
	VectorData_alloc_zero(ptr, n);
	return result;
}

VALUE rb_cVector_of_size(int argc, VALUE *argv, VALUE class) {
	// Parse arguments
	VALUE v_size, v_fill;
	rb_scan_args(argc, argv, "11", &v_size, &v_fill);
	int size = NUM2INT(v_size);
	if (size <= 0) rb_raise(rb_eTypeError, "Vector can not have size 0 or smaller");
	double fill;
	if (v_fill == Qnil) fill = 0.0;	else fill = NUM2DBL(v_fill);
	// Allocate result
	VALUE result = rb_cVector_new(class);
	VectorData *ptr;
	Data_Get_Struct(result,VectorData, ptr);
	VectorData_alloc(ptr, size);
	int i;
	for(i = 0; i < size; i++) {
		ptr->data[i] = fill;
	}
	return result;
}

VALUE rb_cVector_from_array(VALUE class, VALUE array) { ///
	VALUE result = rb_cVector_new(class);
	VectorData *ptr;
	Data_Get_Struct(result,VectorData, ptr);
	Check_Type(array, T_ARRAY); // Check if array is Array
	int size = RARRAY_LEN(array);
	if (size == 0) rb_raise(rb_eTypeError, "Vector can not be created from an empty array");
	VectorData_alloc(ptr, size);

	VALUE* array_ptr = RARRAY_PTR(array);
	int i;
	for(i = 0; i < size; i++) {
		ptr->data[i] = NUM2DBL(array_ptr[i]);
	}

	return result;
}

//==============================================================================
// Vector information
//==============================================================================

VALUE rb_cVector_size(VALUE self) {
	VectorData *ptr;
	Data_Get_Struct(self,VectorData, ptr);
	return INT2NUM(ptr->size);
}

//==============================================================================
// Elemet access
//==============================================================================

VALUE rb_cVector_get_i(VALUE self, VALUE i) {
	VectorData *ptr;
	Data_Get_Struct(self,VectorData, ptr);
	int ii = NUM2INT(i);
	// Out of bounds check
	vector_raise_index_error(ptr, ii);
	return rb_float_new(ptr->data[ii]);
}

VALUE rb_cVector_set_i(VALUE self, VALUE i, VALUE val) {
	VectorData *ptr;
	Data_Get_Struct(self,VectorData, ptr);
	int ii = NUM2INT(i);
	// Out of bounds check
	vector_raise_index_error(ptr, ii);
	ptr->data[ii] = NUM2DBL(val);
	return val;
}

VALUE rb_cVector_subvector(VALUE self, VALUE starti, VALUE size) {
	// Get vector data
	VectorData *ptr;
	Data_Get_Struct(self,VectorData, ptr);
	// Start index and size
	int start_i = NUM2INT(starti);
	int size_i = NUM2INT(size);
	// Error check
	if (size_i <= 0) rb_raise(rb_eTypeError, "Vector can not have size <= 0");
	vector_raise_index_error(ptr, start_i);
	vector_raise_index_error(ptr, start_i + size_i);
	// Allocate result
	VALUE result = rb_cVector_new(rb_class_of(self));
	VectorData *ptr_r;
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, size_i);
	// Memcpy
	memcpy(ptr_r->data, &(ptr->data[start_i]), size_i * sizeof(dbl));
	return result;
}


//==============================================================================
// Misc
//==============================================================================

VALUE rb_cVector_to_a(VALUE self) {
	VectorData *ptr;
	Data_Get_Struct(self,VectorData, ptr);
	VALUE array = rb_ary_new2(ptr->size);
	int i;
	for (i = 0; i < ptr->size; i++) {
		rb_ary_store(array, i, rb_float_new(ptr->data[i]));
	}
	return array;
}

VALUE rb_cVector_to_matrix(VALUE self){
	VectorData *ptr;
	Data_Get_Struct(self,VectorData, ptr);
	// Alocate result
	VALUE result = rb_cMatrix_new(rb_cMatrix);
	MatrixData *ptr_r;
	Data_Get_Struct(result,MatrixData, ptr_r);
	MatrixData_alloc(ptr_r, ptr->size, 1);
	// Memcpy
	memcpy(ptr_r->data, ptr->data, ptr->size * sizeof(dbl));
	return result;
}

VALUE rb_cVector_deep_copy(VALUE self){
	VectorData *ptr;
	Data_Get_Struct(self,VectorData, ptr);
	// Alocate result
	VALUE result = rb_cVector_new(rb_class_of(self));
	VectorData *ptr_r;
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, ptr->size);
	// Memcpy
	memcpy(ptr_r->data, ptr->data, ptr->size * sizeof(dbl));
	return result;
}

//==============================================================================
// Operators
//==============================================================================

VALUE rb_cVector_plus(VALUE self, VALUE vector) {
	vector_raise_if_not_vector(vector);
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	VectorData *ptr_b;
	Data_Get_Struct(vector,VectorData, ptr_b);
	// Check dimensions
	vector_raise_if_not_same_size(ptr_a, ptr_b);
	// Allocate result
	VALUE result = rb_cVector_new(rb_class_of(self));
	VectorData *ptr_r;
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, ptr_a->size);
	// Sum
	int size = ptr_a->size;
	int i;
	for (i = 0; i < size; i++) {
		ptr_r->data[i] = ptr_a->data[i] + ptr_b->data[i];
	}
	// Return
	return result;
}

VALUE rb_cVector_minus(VALUE self, VALUE vector) {
	vector_raise_if_not_vector(vector);
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	VectorData *ptr_b;
	Data_Get_Struct(vector,VectorData, ptr_b);
	// Check dimensions
	vector_raise_if_not_same_size(ptr_a, ptr_b);
	// Allocate result
	VALUE result = rb_cVector_new(rb_class_of(self));
	VectorData *ptr_r;
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, ptr_a->size);
	// Sum
	int size = ptr_a->size;
	int i;
	for (i = 0; i < size; i++) {
		ptr_r->data[i] = ptr_a->data[i] - ptr_b->data[i];
	}
	// Return
	return result;
}

VALUE rb_cVector_plus_inplace(VALUE self, VALUE vector) {
	vector_raise_if_not_vector(vector);
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	VectorData *ptr_b;
	Data_Get_Struct(vector,VectorData, ptr_b);
	// Check dimensions
	vector_raise_if_not_same_size(ptr_a, ptr_b);
	// Sum
	int size = ptr_a->size;
	int i;
	for (i = 0; i < size; i++) {
		ptr_a->data[i] = ptr_a->data[i] + ptr_b->data[i];
	}
	// Return
	return Qnil;
}

VALUE rb_cVector_minus_inplace(VALUE self, VALUE vector) {
	vector_raise_if_not_vector(vector);
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	VectorData *ptr_b;
	Data_Get_Struct(vector,VectorData, ptr_b);
	// Check dimensions
	vector_raise_if_not_same_size(ptr_a, ptr_b);
	// Sum
	int size = ptr_a->size;
	int i;
	for (i = 0; i < size; i++) {
		ptr_a->data[i] = ptr_a->data[i] - ptr_b->data[i];
	}
	// Return
	return Qnil;
}

VALUE rb_cVector_comutative_multiply(VALUE self, VALUE number) {
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	double n = NUM2DBL(number);
	// Allocate result
	VALUE result = rb_cVector_new(rb_class_of(self));
	VectorData *ptr_r;
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, ptr_a->size);
	// Traverse vector
	int size = ptr_a->size;
	int i;
	for (i = 0; i < size; i++) {
		ptr_r->data[i] = ptr_a->data[i] * n;
	}
	// Return
	return result;
}

VALUE rb_cVector_multiply_inplace(VALUE self, VALUE number) {
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	double n = NUM2DBL(number);
	// Traverse vector
	int size = ptr_a->size;
	int i;
	for (i = 0; i < size; i++) {
		ptr_a->data[i] = ptr_a->data[i] * n;
	}
	// Return
	return Qnil;
}

VALUE rb_cVector_dot(VALUE self, VALUE vector) {
	vector_raise_if_not_vector(vector);
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	VectorData *ptr_b;
	Data_Get_Struct(vector,VectorData, ptr_b);
	// Check dimensions
	vector_raise_if_not_same_size(ptr_a, ptr_b);
	// Sum
	int size = ptr_a->size;
	int i;
	double sum = 0.0;
	for (i = 0; i < size; i++) {
		sum += ptr_a->data[i] * ptr_b->data[i];
	}
	// Return
	return rb_float_new(sum);
}

VALUE rb_cVector_elementwise_multiply(VALUE self, VALUE vector) {
	vector_raise_if_not_vector(vector);
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	VectorData *ptr_b;
	Data_Get_Struct(vector,VectorData, ptr_b);
	// Check dimensions
	vector_raise_if_not_same_size(ptr_a, ptr_b);
	// Allocate result
	VALUE result = rb_cVector_new(rb_class_of(self));
	VectorData *ptr_r;
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, ptr_a->size);
	// Sum
	int size = ptr_a->size;
	int i;
	for (i = 0; i < size; i++) {
		ptr_r->data[i] = ptr_a->data[i] * ptr_b->data[i];
	}
	// Return
	return result;
}

VALUE rb_cVector_elementwise_divide(VALUE self, VALUE vector) {
	vector_raise_if_not_vector(vector);
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	VectorData *ptr_b;
	Data_Get_Struct(vector,VectorData, ptr_b);
	// Check dimensions
	vector_raise_if_not_same_size(ptr_a, ptr_b);
	// Allocate result
	VALUE result = rb_cVector_new(rb_class_of(self));
	VectorData *ptr_r;
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, ptr_a->size);
	// Sum
	int size = ptr_a->size;
	int i;
	for (i = 0; i < size; i++) {
		ptr_r->data[i] = ptr_a->data[i] / ptr_b->data[i];
	}
	// Return
	return result;
}

//==============================================================================
// Comparison
//==============================================================================

VALUE rb_cVector_equal(VALUE self, VALUE vector) {
	if (CLASS_OF(vector) != CLASS_OF(self)) return Qfalse;
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	VectorData *ptr_b;
	Data_Get_Struct(vector,VectorData, ptr_b);
	// Check dimensions
	vector_raise_if_not_same_size(ptr_a, ptr_b);
	// Sum
	int size = ptr_a->size;
	int i;
	for (i = 0; i < size; i++) {
		if (ptr_a->data[i] != ptr_b->data[i]) return Qfalse;
	}
	return Qtrue;
}

VALUE rb_cVector_equal_eps(VALUE self, VALUE vector) {
	if (CLASS_OF(vector) != CLASS_OF(self)) return Qfalse;
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	VectorData *ptr_b;
	Data_Get_Struct(vector,VectorData, ptr_b);
	// Check dimensions
	vector_raise_if_not_same_size(ptr_a, ptr_b);
	// Sum
	int size = ptr_a->size;
	int i;
	dbl epsilon = NUM2DBL(rb_cv_get(rb_cVector, "@@default_epsilon"));
	for (i = 0; i < size; i++) {
		if (fabs(ptr_a->data[i] - ptr_b->data[i]) > epsilon) return Qfalse;
	}
	return Qtrue;
}

//==============================================================================
// Vector operations
//==============================================================================

VALUE rb_cVector_sum(VALUE self) {
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	// Sum
	int size = ptr_a->size;
	int i;
	double sum = 0.0;
	for (i = 0; i < size; i++) {
		sum += ptr_a->data[i];
	}
	// Return
	return rb_float_new(sum);
}

VALUE rb_cVector_abs(VALUE self) {
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	// Sum
	int size = ptr_a->size;
	int i;
	double sum = 0.0;
	for (i = 0; i < size; i++) {
		sum += ptr_a->data[i] * ptr_a->data[i];
	}
	// Return
	return rb_float_new(sqrt(sum));
}

VALUE rb_cVector_max(VALUE self) {
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	// Sum
	int size = ptr_a->size;
	int i;
	double m = ptr_a->data[0];
	for (i = 0; i < size; i++) {
		if (ptr_a->data[i] > m) m = ptr_a->data[i];
	}
	// Return
	return rb_float_new(m);
}

VALUE rb_cVector_min(VALUE self) {
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	// Sum
	int size = ptr_a->size;
	int i;
	double m = ptr_a->data[0];
	for (i = 0; i < size; i++) {
		if (ptr_a->data[i] < m) m = ptr_a->data[i];
	}
	// Return
	return rb_float_new(m);
}

VALUE rb_cVector_max_abs(VALUE self) {
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	// Sum
	int size = ptr_a->size;
	int i;
	double m = fabs(ptr_a->data[0]);
	double a;
	for (i = 0; i < size; i++) {
		a = fabs(ptr_a->data[i]);
		if (a > m) m = a;
	}
	// Return
	return rb_float_new(m);
}

VALUE rb_cVector_min_abs(VALUE self) {
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	// Sum
	int size = ptr_a->size;
	int i;
	double m = fabs(ptr_a->data[0]);
	double a;
	for (i = 0; i < size; i++) {
		a = fabs(ptr_a->data[i]);
		if (a < m) m = a;
	}
	// Return
	return rb_float_new(m);
}

VALUE rb_cVector_normalize(VALUE self) {
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	// Allocate result
	VALUE result = rb_cVector_new(rb_class_of(self));
	VectorData *ptr_r;
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, ptr_a->size);
	// get size
	int size = ptr_a->size;
	int i;
	double s = 0.0;
	for (i = 0; i < size; i++) {
		s = s + ptr_a->data[i] * ptr_a->data[i];
	}
	s = sqrt(s);
	// build result
	for (i = 0; i < size; i++) {
		ptr_r->data[i] = ptr_a->data[i] / s;
	}
	// Return
	return result;
}

VALUE rb_cVector_normalize_inplace(VALUE self) {
	VectorData *ptr_a;
	Data_Get_Struct(self,VectorData, ptr_a);
	// get size
	int size = ptr_a->size;
	int i;
	double s = 0.0;
	for (i = 0; i < size; i++) {
		s += ptr_a->data[i] * ptr_a->data[i];
	}
	s = sqrt(s);
	// resize
	for (i = 0; i < size; i++) {
		ptr_a->data[i] = ptr_a->data[i] / s;
	}
	// Return
	return Qnil;
}

//==============================================================================
// Marshalling
//==============================================================================

VALUE rb_cVector_marshal_dump(VALUE self, VALUE depth) {
	// Get vector data
	VectorData *ptr;
	Data_Get_Struct(self,VectorData, ptr);
	// Create string
	VALUE str = rb_str_buf_new(sizeof(int) + ptr->size * sizeof(dbl));
	// Wrire string
	rb_str_buf_cat(str, (const char*)(&(ptr->size)), sizeof(int)) ;
	rb_str_buf_cat(str, (const char*)(ptr->data), ptr->size * sizeof(dbl));
	return str;
}

VALUE rb_cVector_marshal_load(VALUE class, VALUE str) {
	// Convert value to C string
	char* p = (char*)RSTRING_PTR(str);
	// Get matrix size
	int* n = (int*)p;
	// Allocate result
	VALUE result = rb_cVector_new(class);
	VectorData *ptr_r;
	Data_Get_Struct(result,VectorData, ptr_r);
	VectorData_alloc(ptr_r, n[0]);
	// Copy data
	memcpy(ptr_r->data, p + sizeof(int), ptr_r->size * sizeof(dbl));
	return result;
}

//##############################################################################
// Extension init
//##############################################################################

void Init_vector_c () {

	rb_cVector = rb_define_class_under(rb_mAlgebra, "Vector", rb_cObject);

	rb_define_alloc_func(rb_cVector, rb_cVector_new);

	// Constructors
	rb_define_method(rb_cVector, "initialize", rb_cVector_initialize, 1);
	rb_define_singleton_method(rb_cVector, "from_array", rb_cVector_from_array, 1);
	rb_define_singleton_method(rb_cVector, "of_size", rb_cVector_of_size, -1);
	rb_define_singleton_method(rb_cVector, "zero", rb_cVector_zero, 1);

	// Vector information
	rb_define_method(rb_cVector,"size", rb_cVector_size, 0);

	// Elemet access
	rb_define_method(rb_cVector,"[]", rb_cVector_get_i, 1);
	rb_define_method(rb_cVector,"[]=", rb_cVector_set_i, 2);
	rb_define_method(rb_cVector,"subvector", rb_cVector_subvector, 2);

	// Misc
	rb_define_method(rb_cVector, "to_a", rb_cVector_to_a, 0);
	rb_define_method(rb_cVector, "to_matrix", rb_cVector_to_matrix, 0);
	rb_define_method(rb_cVector, "deep_copy", rb_cVector_deep_copy, 0);

	// Operators
	rb_define_method(rb_cVector, "+", rb_cVector_plus, 1);
	rb_define_method(rb_cVector, "-", rb_cVector_minus, 1);
	rb_define_method(rb_cVector, "plus!", rb_cVector_plus_inplace, 1);
	rb_define_method(rb_cVector, "minus!", rb_cVector_minus_inplace, 1);
	rb_define_method(rb_cVector, "comutative_multiply", rb_cVector_comutative_multiply, 1);
	rb_define_method(rb_cVector, "mul!", rb_cVector_multiply_inplace, 1);
	rb_define_method(rb_cVector, "dot", rb_cVector_dot, 1);
	rb_define_method(rb_cVector, "elementwise_multiply", rb_cVector_elementwise_multiply, 1);
	rb_define_method(rb_cVector, "elementwise_divide", rb_cVector_elementwise_divide, 1);

	// Comparison
	rb_define_method(rb_cVector, "==", rb_cVector_equal, 1);
	rb_define_method(rb_cVector, "=~", rb_cVector_equal_eps, 1);

	// Vector operations
	rb_define_method(rb_cVector, "sum", rb_cVector_sum, 0);
	rb_define_method(rb_cVector, "abs", rb_cVector_abs, 0);
	rb_define_method(rb_cVector, "max", rb_cVector_max, 0);
	rb_define_method(rb_cVector, "min", rb_cVector_min, 0);
	rb_define_method(rb_cVector, "max_abs", rb_cVector_max_abs, 0);
	rb_define_method(rb_cVector, "min_abs", rb_cVector_min_abs, 0);
	rb_define_method(rb_cVector, "normalize", rb_cVector_normalize, 0);
	rb_define_method(rb_cVector, "normalize!", rb_cVector_normalize_inplace, 0);

	// Marshaling
	rb_define_method(rb_cVector, "_dump", RUBY_METHOD_FUNC(rb_cVector_marshal_dump), 1);
	rb_define_singleton_method(rb_cVector, "_load", RUBY_METHOD_FUNC(rb_cVector_marshal_load), 1);

	///### Check type where argument is Vector
	///### Check for zero size at init

}
