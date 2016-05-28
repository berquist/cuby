#include "ruby.h"
#include "math.h"


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


static VALUE rb_cCoordinate;

typedef struct coordinatedata {
	double x;
	double y;
	double z;
} CoordinateData;

CoordinateData *CoordinateDataNew() {
	CoordinateData *c = ALLOC(CoordinateData);
	//CoordinateData *c = malloc(sizeof(CoordinateData));
	return c;
}

CoordinateData *CoordinateDataNewZero() {
	CoordinateData *c = ALLOC(CoordinateData);
	//CoordinateData *c = malloc(sizeof(CoordinateData));
	c->x = 0.0;
	c->y = 0.0;
	c->z = 0.0;
	return c;
}

CoordinateData *CoordinateDataNewXYZ(double x, double y, double z) {
	CoordinateData *c = ALLOC(CoordinateData);
	//CoordinateData *c = malloc(sizeof(CoordinateData));
	c->x = x;
	c->y = y;
	c->z = z;
	return c;
}

void CoordinateDataFree(CoordinateData *ptr){
	//printf("Ptr: %X\n",(int)ptr);
	free(ptr);
}

//------------------------------------------------------------------------------

static VALUE coord_new_noinit(VALUE class) {
	CoordinateData *ptr = CoordinateDataNew();
	VALUE val = Data_Wrap_Struct(class, 0, CoordinateDataFree, ptr);
	return val;
}

static VALUE coord_initialize(int argc, VALUE *argv, VALUE self) {
	double x, y, z;
	CoordinateData *p = DATA_PTR(self);
	//CoordinateData *p;
	//Data_Get_Struct(self,CoordinateData, p);
	if (argc > 0) {x = NUM2DBL(argv[0]);} else {x = 0.0;}
	if (argc > 1) {y = NUM2DBL(argv[1]);} else {y = 0.0;}
	if (argc > 2) {z = NUM2DBL(argv[2]);} else {z = 0.0;}
	
	p->x = x;
	p->y = y;
	p->z = z;
	return self;
}

static VALUE coord_from_values(int argc, VALUE *argv, VALUE class) {
	double x, y, z;
	if (argc > 0) {x = NUM2DBL(argv[0]);} else {x = 0.0;}
	if (argc > 1) {y = NUM2DBL(argv[1]);} else {y = 0.0;}
	if (argc > 2) {z = NUM2DBL(argv[2]);} else {z = 0.0;}
	
	CoordinateData *ptr = CoordinateDataNewXYZ(x,y,z);
	VALUE data = Data_Wrap_Struct(class, 0, CoordinateDataFree, ptr);
	return data;
}

static VALUE coord_zero(VALUE class) {
	CoordinateData *ptr = CoordinateDataNewZero();
	VALUE data = Data_Wrap_Struct(class, 0, CoordinateDataFree, ptr);
	return data;
}

static VALUE coord_from_array(VALUE class, VALUE array) {
	VALUE* a = RARRAY_PTR(array);
	CoordinateData *ptr = CoordinateDataNewXYZ(NUM2DBL(a[0]), NUM2DBL(a[1]), NUM2DBL(a[2]));
	VALUE data = Data_Wrap_Struct(class, 0, CoordinateDataFree, ptr);
	return data;
}

static VALUE coord_from_coord(VALUE class, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Coordinate can be created only from a Coordinate-like type");
	CoordinateData *o;
	Data_Get_Struct(other,CoordinateData, o);
	CoordinateData *ptr = CoordinateDataNewXYZ(o->x,o->y,o->z);
	VALUE data = Data_Wrap_Struct(class, 0, CoordinateDataFree, ptr);
	return data;
}

//------------------------------------------------------------------------------
// Access to elements

static VALUE coord_x(VALUE self) {
	CoordinateData *p;
	Data_Get_Struct(self,CoordinateData, p);
	return rb_float_new(p->x);
}

static VALUE coord_y(VALUE self) {
	CoordinateData *p;
	Data_Get_Struct(self,CoordinateData, p);
	return rb_float_new(p->y);
}

static VALUE coord_z(VALUE self) {
	CoordinateData *p;
	Data_Get_Struct(self,CoordinateData, p);
	return rb_float_new(p->z);
}

static VALUE coord_x_assign(VALUE self, VALUE val) {
	CoordinateData *p;
	Data_Get_Struct(self,CoordinateData, p);
	p->x = NUM2DBL(val);
	return rb_float_new(p->x);
}

static VALUE coord_y_assign(VALUE self, VALUE val) {
	CoordinateData *p;
	Data_Get_Struct(self,CoordinateData, p);
	p->y = NUM2DBL(val);
	return rb_float_new(p->y);
}

static VALUE coord_z_assign(VALUE self, VALUE val) {
	CoordinateData *p;
	Data_Get_Struct(self,CoordinateData, p);
	p->z = NUM2DBL(val);
	return rb_float_new(p->z);
}

//------------------------------------------------------------------------------

static VALUE coord_set_coord(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Coordinates can be set ofnly from a Coordinate-like type");
	CoordinateData *s, *o;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	s->x = o->x;
	s->y = o->y;
	s->z = o->z;
	return Qnil;
}

//------------------------------------------------------------------------------
// Array-like access

static VALUE coord_at_index(VALUE self, VALUE index) {
	CoordinateData *s;
	int i = FIX2INT(index);
	// check index range
	if (i > 2 || i < 0) rb_raise(rb_eIndexError, "Index not in range 0..2");
	Data_Get_Struct(self,CoordinateData, s);
	return rb_float_new(((double*)s)[i]);
}

static VALUE coord_assign_at_index(VALUE self, VALUE index, VALUE val) {
	CoordinateData *s;
	int i = FIX2INT(index);
	// check index range
	if (i > 2 || i < 0) rb_raise(rb_eIndexError, "Index not in range 0..2");
	Data_Get_Struct(self,CoordinateData, s);
	((double*)s)[i] = NUM2DBL(val);
	return rb_float_new(((double*)s)[i]);
}

//------------------------------------------------------------------------------
// Type conversion

static VALUE coord_to_s(VALUE self) {
	CoordinateData *s;
	char string[100];
	Data_Get_Struct(self,CoordinateData, s);
	sprintf(string, "[%f, %f, %f]",s->x,s->y,s->z);
	return rb_str_new2(string);
}

static VALUE coord_inspect(VALUE self) {
	CoordinateData *s;
	char string[100];
	Data_Get_Struct(self,CoordinateData, s);
	sprintf(string, "\n[%f, %f, %f]",s->x,s->y,s->z);
	return rb_str_new2(string);
}

static VALUE coord_to_a(VALUE self) {
	CoordinateData *s;
	Data_Get_Struct(self,CoordinateData, s);
	return rb_ary_new3(3, rb_float_new(s->x), rb_float_new(s->y), rb_float_new(s->z)); 
}

static VALUE coord_to_coordinate(VALUE self) {
	CoordinateData *s;
	Data_Get_Struct(self,CoordinateData, s);
	CoordinateData *ptr = CoordinateDataNewXYZ(s->x,s->y,s->z);
	VALUE data = Data_Wrap_Struct(rb_cCoordinate, 0, CoordinateDataFree, ptr);
	return data;
}

//------------------------------------------------------------------------------
// Iterators

static VALUE coord_each_index(VALUE self) {
	rb_yield(INT2FIX(0));
	rb_yield(INT2FIX(1));
	rb_yield(INT2FIX(2));
	return Qnil;
}

static VALUE coord_each(VALUE self) {
	CoordinateData *s;
	Data_Get_Struct(self,CoordinateData, s);
	rb_yield(rb_float_new(s->x));
	rb_yield(rb_float_new(s->y));
	rb_yield(rb_float_new(s->z));
	return Qnil;
}

static VALUE coord_map(VALUE self) {
	CoordinateData *s;
	Data_Get_Struct(self,CoordinateData, s);
	double xx = NUM2DBL(rb_yield(rb_float_new(s->x)));
	double yy = NUM2DBL(rb_yield(rb_float_new(s->y)));
	double zz = NUM2DBL(rb_yield(rb_float_new(s->z)));
	CoordinateData *ptr = CoordinateDataNewXYZ(xx,yy,zz);
	VALUE data = Data_Wrap_Struct(rb_cCoordinate, 0, CoordinateDataFree, ptr);
	return data;
}

static VALUE coord_map_inplace(VALUE self) {
	CoordinateData *s;
	Data_Get_Struct(self,CoordinateData, s);
	s->x = NUM2DBL(rb_yield(rb_float_new(s->x)));
	s->y = NUM2DBL(rb_yield(rb_float_new(s->y)));
	s->z = NUM2DBL(rb_yield(rb_float_new(s->z)));
	return Qnil;
}

//------------------------------------------------------------------------------

static VALUE coord_unary_plus(VALUE self) {
	return self;
}

static VALUE coord_unary_minus(VALUE self) {
	CoordinateData *s, *r;
	Data_Get_Struct(self,CoordinateData, s);
	r = CoordinateDataNew();
	r->x = s->x * -1.0;
	r->y = s->y * -1.0;
	r->z = s->z * -1.0;
	VALUE result = Data_Wrap_Struct(rb_cCoordinate, 0, CoordinateDataFree, r);
	return result;
}

//------------------------------------------------------------------------------

static VALUE coord_plus_inplace(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Operand must be a Coordinate-like type");
	CoordinateData *s, *o;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	s->x += o->x;
	s->y += o->y;
	s->z += o->z;
	return self;
}

static VALUE coord_minus_inplace(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Operand must be a Coordinate-like type");
	CoordinateData *s, *o;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	s->x -= o->x;
	s->y -= o->y;
	s->z -= o->z;
	return self;
}

//------------------------------------------------------------------------------

static VALUE coord_plus(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Operand must be a Coordinate-like type");
	CoordinateData *s, *o, *r;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	r = CoordinateDataNew();
	r->x = s->x + o->x;
	r->y = s->y + o->y;
	r->z = s->z + o->z;
	VALUE result = Data_Wrap_Struct(rb_cCoordinate, 0, CoordinateDataFree, r);
	return result;
}

static VALUE coord_minus(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Operand must be a Coordinate-like type");
	CoordinateData *s, *o, *r;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	r = CoordinateDataNew();
	r->x = s->x - o->x;
	r->y = s->y - o->y;
	r->z = s->z - o->z;
	VALUE result = Data_Wrap_Struct(rb_cCoordinate, 0, CoordinateDataFree, r);
	return result;
}

//------------------------------------------------------------------------------

static VALUE coord_dot(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Operand must be a Coordinate-like type");
	CoordinateData *s, *o;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	return rb_float_new(s->x*o->x + s->y*o->y + s->z*o->z);
}

static VALUE coord_cross(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Operand must be a Coordinate-like type");
	CoordinateData *s, *o, *r;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	r = CoordinateDataNew();
	r->x = s->y * o->z - s->z * o->y;
	r->y = s->z * o->x - s->x * o->z;
	r->z = s->x * o->y - s->y * o->x;
	VALUE result = Data_Wrap_Struct(rb_cCoordinate, 0, CoordinateDataFree, r);
	return result;
}

//------------------------------------------------------------------------------

static VALUE coord_mul_number(VALUE self, VALUE other) {
	CoordinateData *s, *r;
	double o = NUM2DBL(other);
	Data_Get_Struct(self,CoordinateData, s);
	r = CoordinateDataNew();
	r->x = s->x * o;
	r->y = s->y * o;
	r->z = s->z * o;
	VALUE result = Data_Wrap_Struct(rb_cCoordinate, 0, CoordinateDataFree, r);
	return result;
}

static VALUE coord_mul(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other,rb_cNumeric) == Qtrue) return coord_mul_number(self,other);
	if (rb_obj_is_kind_of(other,rb_cCoordinate) == Qtrue) return coord_dot(self,other);
	rb_raise(rb_eTypeError, "Coordinate can be multiplied by coordinate or number only");
}

static VALUE coord_div(VALUE self, VALUE other) {
	CoordinateData *s, *r;
	double o = NUM2DBL(other);
	Data_Get_Struct(self,CoordinateData, s);
	r = CoordinateDataNew();
	r->x = s->x / o;
	r->y = s->y / o;
	r->z = s->z / o;
	VALUE result = Data_Wrap_Struct(rb_cCoordinate, 0, CoordinateDataFree, r);
	return result;
}

//------------------------------------------------------------------------------
// Comparison

static VALUE coord_equal(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other,rb_cCoordinate) != Qtrue) return Qfalse;
	CoordinateData *s, *o;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	if (s->x == o->x && s->y == o->y &&  s->z == o->z) return Qtrue;
	return Qfalse;
}

static VALUE coord_similar(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other,rb_cCoordinate) != Qtrue) return Qfalse;
	CoordinateData *s, *o;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	double epsilon = NUM2DBL(rb_cv_get(rb_cCoordinate, "@@epsilon")); 
	if (fabs(s->x - o->x) < epsilon && fabs(s->y - o->y) < epsilon && fabs(s->z - o->z) < epsilon) return Qtrue;
	return Qfalse;
}

static VALUE is_zero(VALUE self) {
	CoordinateData *s;
	Data_Get_Struct(self,CoordinateData, s);
	if (s->x == 0.0 && s->y == 0.0 &&  s->z == 0.0) return Qtrue;
	return Qfalse;
}

static VALUE coord_get_epsilon(VALUE class) {
	return rb_cv_get(rb_cCoordinate, "@@epsilon");
}

static VALUE coord_set_epsilon(VALUE class, VALUE val) {
	rb_cv_set(rb_cCoordinate, "@@epsilon", val);
	return val;
}

//------------------------------------------------------------------------------
// Misc

static VALUE coord_abs(VALUE self) {
	CoordinateData *s;
	Data_Get_Struct(self,CoordinateData, s);
	return  rb_float_new(sqrt(s->x*s->x + s->y*s->y + s->z*s->z));
}

static VALUE coord_normalize(VALUE self) {
	CoordinateData *s, *r;
	Data_Get_Struct(self,CoordinateData, s);
	r = CoordinateDataNew();
	double abs = sqrt(s->x*s->x + s->y*s->y + s->z*s->z);
	r->x = s->x / abs;
	r->y = s->y / abs;
	r->z = s->z / abs;
	VALUE result = Data_Wrap_Struct(rb_cCoordinate, 0, CoordinateDataFree, r);
	return result;
}

static VALUE coord_normalize_inplace(VALUE self) {
	CoordinateData *s;
	Data_Get_Struct(self,CoordinateData, s);
	double abs = sqrt(s->x*s->x + s->y*s->y + s->z*s->z);
	s->x = s->x / abs;
	s->y = s->y / abs;
	s->z = s->z / abs;
	return Qnil;
}

static VALUE coord_distance(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Argument must be a Coordinate-like type");
	CoordinateData *s, *o;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	double dx = s->x - o->x;
	double dy = s->y - o->y;
	double dz = s->z - o->z;
	return rb_float_new(sqrt(dx*dx + dy*dy + dz*dz));
}

static VALUE coord_angle(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Argument must be a Coordinate-like type");
	CoordinateData *s, *o;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	double abs1 = sqrt(s->x*s->x + s->y*s->y + s->z*s->z);
	double abs2 = sqrt(o->x*o->x + o->y*o->y + o->z*o->z);
	if (abs1 == 0.0) rb_raise(rb_eArgError, "Coordinate for angle calculation can not be zero");
	if (abs2 == 0.0) rb_raise(rb_eArgError, "Coordinate for angle calculation can not be zero");
	double dot = s->x*o->x + s->y*o->y + s->z*o->z;
	double cos = dot/abs1/abs2;
	// Numerical issue: can be close to 1 but out of interval:
	if (cos < -1.0) cos = -1.0;
	if (cos > 1.0) cos = 1.0;
	return rb_float_new(acos(cos));
}

static VALUE coord_angle_deg(VALUE self, VALUE other) {
	if (rb_obj_is_kind_of(other, rb_cCoordinate) != Qtrue) rb_raise(rb_eTypeError, "Argument must be a Coordinate-like type");
	CoordinateData *s, *o;
	Data_Get_Struct(self,CoordinateData, s);
	Data_Get_Struct(other,CoordinateData, o);
	double abs1 = sqrt(s->x*s->x + s->y*s->y + s->z*s->z);
	double abs2 = sqrt(o->x*o->x + o->y*o->y + o->z*o->z);
	double dot = s->x*o->x + s->y*o->y + s->z*o->z;
	return rb_float_new(acos(dot/abs1/abs2) * 180.0 / M_PI);
}

static VALUE coord_marshal_dump(VALUE self, VALUE depth) {
	CoordinateData *s;
	Data_Get_Struct(self,CoordinateData, s);
	VALUE str = rb_str_buf_new(1*sizeof(CoordinateData));
	rb_str_buf_cat(str, (const char*)(s), 1*sizeof(CoordinateData)) ;
	return str;
}

static VALUE coord_marshal_load(VALUE class, VALUE str) {
	//char* p = (char*)RSTRING(str)->ptr;
	char* p = (char*)RSTRING_PTR(str);
	CoordinateData *s = (CoordinateData*)p;
	CoordinateData *ptr = CoordinateDataNewXYZ(s->x,s->y,s->z);
	VALUE data = Data_Wrap_Struct(class, 0, CoordinateDataFree, ptr);
	return data;
}

//------------------------------------------------------------------------------
// Deep copy

static VALUE coord_deep_copy(VALUE self) {
	CoordinateData *s;
	Data_Get_Struct(self,CoordinateData, s);
	CoordinateData *ptr = CoordinateDataNewXYZ(s->x,s->y,s->z);
	VALUE data = Data_Wrap_Struct(rb_cCoordinate, 0, CoordinateDataFree, ptr);
	return data;
}

//==============================================================================
// Version information (returns :c)
//==============================================================================
VALUE rb_cCoordinate_version(VALUE self) {
	return ID2SYM(rb_intern("c"));
}

/*------------------------------------------------------------------------------
 * Coordinate class
 *----------------------------------------------------------------------------*/


void Init_coordinate_c() {
	rb_cCoordinate = rb_define_class("Coordinate", rb_cObject);

	rb_define_singleton_method(rb_cCoordinate, "version", rb_cCoordinate_version, 0);

	// Class variables
	rb_define_class_variable(rb_cCoordinate, "@@epsilon", rb_float_new(1.0e-8));

	// Alocation function
	rb_define_alloc_func(rb_cCoordinate, coord_new_noinit);

	// Constructors
	rb_define_method(rb_cCoordinate, "initialize", RUBY_METHOD_FUNC(coord_initialize), -1);
	rb_define_module_function(rb_cCoordinate, "[]", RUBY_METHOD_FUNC(coord_from_values), -1);
	rb_define_module_function(rb_cCoordinate, "zero", RUBY_METHOD_FUNC(coord_zero), 0);
	rb_define_module_function(rb_cCoordinate, "from_array", RUBY_METHOD_FUNC(coord_from_array), 1);
	rb_define_module_function(rb_cCoordinate, "from_coordinate", RUBY_METHOD_FUNC(coord_from_coord), 1);
	
	// Element access
	rb_define_method(rb_cCoordinate, "x", RUBY_METHOD_FUNC(coord_x), 0);
	rb_define_method(rb_cCoordinate, "y", RUBY_METHOD_FUNC(coord_y), 0);
	rb_define_method(rb_cCoordinate, "z", RUBY_METHOD_FUNC(coord_z), 0);
	rb_define_method(rb_cCoordinate, "x=", RUBY_METHOD_FUNC(coord_x_assign), 1);
	rb_define_method(rb_cCoordinate, "y=", RUBY_METHOD_FUNC(coord_y_assign), 1);
	rb_define_method(rb_cCoordinate, "z=", RUBY_METHOD_FUNC(coord_z_assign), 1);
	rb_define_method(rb_cCoordinate, "set_coord", RUBY_METHOD_FUNC(coord_set_coord), 1);
	rb_define_method(rb_cCoordinate, "[]", RUBY_METHOD_FUNC(coord_at_index), 1);
	rb_define_method(rb_cCoordinate, "[]=", RUBY_METHOD_FUNC(coord_assign_at_index), 2);

	// Type conversion
	rb_define_method(rb_cCoordinate, "to_s", RUBY_METHOD_FUNC(coord_to_s), 0);
	rb_define_method(rb_cCoordinate, "inspect", RUBY_METHOD_FUNC(coord_inspect), 0);
	rb_define_method(rb_cCoordinate, "to_a", RUBY_METHOD_FUNC(coord_to_a), 0);
	rb_define_method(rb_cCoordinate, "to_coordinate", RUBY_METHOD_FUNC(coord_to_coordinate), 0);

	// Iterators
	rb_define_method(rb_cCoordinate, "each_index", RUBY_METHOD_FUNC(coord_each_index), 0);
	rb_define_method(rb_cCoordinate, "each", RUBY_METHOD_FUNC(coord_each), 0);
	rb_define_method(rb_cCoordinate, "map", RUBY_METHOD_FUNC(coord_map), 0);
	rb_define_method(rb_cCoordinate, "collect", RUBY_METHOD_FUNC(coord_map), 0);
	rb_define_method(rb_cCoordinate, "map!", RUBY_METHOD_FUNC(coord_map_inplace), 0);

	// Operators
	rb_define_method(rb_cCoordinate, "+@", RUBY_METHOD_FUNC(coord_unary_plus), 0);
	rb_define_method(rb_cCoordinate, "-@", RUBY_METHOD_FUNC(coord_unary_minus), 0);
	rb_define_method(rb_cCoordinate, "plus!", RUBY_METHOD_FUNC(coord_plus_inplace), 1);
	rb_define_method(rb_cCoordinate, "minus!", RUBY_METHOD_FUNC(coord_minus_inplace), 1);
	rb_define_method(rb_cCoordinate, "+", RUBY_METHOD_FUNC(coord_plus), 1);
	rb_define_method(rb_cCoordinate, "-", RUBY_METHOD_FUNC(coord_minus), 1);
	rb_define_method(rb_cCoordinate, "dot", RUBY_METHOD_FUNC(coord_dot), 1);
	rb_define_method(rb_cCoordinate, "dot_product", RUBY_METHOD_FUNC(coord_dot), 1);
	rb_define_method(rb_cCoordinate, "cross_product", RUBY_METHOD_FUNC(coord_cross), 1);
	rb_define_method(rb_cCoordinate, "*", RUBY_METHOD_FUNC(coord_mul), 1);
	rb_define_method(rb_cCoordinate, "/", RUBY_METHOD_FUNC(coord_div), 1);

	// Comparison
	rb_define_method(rb_cCoordinate, "==", RUBY_METHOD_FUNC(coord_equal), 1);
	rb_define_method(rb_cCoordinate, "=~", RUBY_METHOD_FUNC(coord_similar), 1);
	rb_define_method(rb_cCoordinate, "zero?", RUBY_METHOD_FUNC(is_zero), 0);
	rb_define_module_function(rb_cCoordinate, "epsilon", RUBY_METHOD_FUNC(coord_get_epsilon), 0);
	rb_define_module_function(rb_cCoordinate, "epsilon=", RUBY_METHOD_FUNC(coord_set_epsilon), 1);
	
	// Misc
	rb_define_method(rb_cCoordinate, "abs", RUBY_METHOD_FUNC(coord_abs), 0);
	rb_define_method(rb_cCoordinate, "absolute", RUBY_METHOD_FUNC(coord_abs), 0);
	rb_define_method(rb_cCoordinate, "r", RUBY_METHOD_FUNC(coord_abs), 0);
	rb_define_method(rb_cCoordinate, "normalize", RUBY_METHOD_FUNC(coord_normalize), 0);
	rb_define_method(rb_cCoordinate, "normalize!", RUBY_METHOD_FUNC(coord_normalize_inplace), 0);
	rb_define_method(rb_cCoordinate, "distance", RUBY_METHOD_FUNC(coord_distance), 1);
	rb_define_method(rb_cCoordinate, "dist", RUBY_METHOD_FUNC(coord_distance), 1);
	rb_define_method(rb_cCoordinate, "angle", RUBY_METHOD_FUNC(coord_angle), 1);
	rb_define_method(rb_cCoordinate, "angle_deg", RUBY_METHOD_FUNC(coord_angle_deg), 1);

	// Marshaling
	rb_define_method(rb_cCoordinate, "_dump", RUBY_METHOD_FUNC(coord_marshal_dump), 1);
	rb_define_singleton_method(rb_cCoordinate, "_load", RUBY_METHOD_FUNC(coord_marshal_load), 1);

	// Deep copy
	rb_define_method(rb_cCoordinate, "deep_copy", RUBY_METHOD_FUNC(coord_deep_copy), 0);
	rb_define_method(rb_cCoordinate, "dup", RUBY_METHOD_FUNC(coord_deep_copy), 0);
}
