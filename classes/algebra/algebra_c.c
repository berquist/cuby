#include <ruby.h>
#include <math.h>

#include "algebra_c.h"

VALUE rb_eDimensionError;
VALUE rb_eSparseWriteError;
VALUE rb_eSparseTypeError;
VALUE rb_mAlgebra;

//==============================================================================
// Version information (returns :c)
//==============================================================================
VALUE rb_mAlgebra_version(VALUE self) {
	return ID2SYM(rb_intern("c"));
}


//==============================================================================
// Module initialization
//==============================================================================
void Init_algebra_c() {
	// Algebra module
	rb_mAlgebra = rb_define_module("Algebra");
	rb_define_singleton_method(rb_mAlgebra,"version", rb_mAlgebra_version, 0);


	// Exceptions
	rb_eDimensionError = rb_define_class_under(rb_mAlgebra, "DimensionError", rb_eRuntimeError);
	rb_eSparseWriteError = rb_define_class_under(rb_mAlgebra, "SparseWriteError", rb_eRuntimeError);
	rb_eSparseTypeError = rb_define_class_under(rb_mAlgebra, "SparseTypeError", rb_eRuntimeError);

	// Initialize classes
	Init_vector_c();
	Init_matrix_c();
	Init_sparse_matrix_csx_c();
}
