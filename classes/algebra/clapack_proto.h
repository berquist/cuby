// Adapted from http://www.netlib.org/clapack/clapack.h

#ifndef __CLAPACK_PROTO_H
#define __CLAPACK_PROTO_H

extern int dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, 
	double *b, int *ldb, double *beta, double *c__, int *ldc);

extern int dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, 
	int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);

extern int dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);

extern int dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u,
	int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);

extern int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

extern int dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

extern int dpotrf_(char *uplo, int *n, double *a, int *lda, int *info); // Cholesky decomposition


#endif /* __CLAPACK_PROTO_H */
