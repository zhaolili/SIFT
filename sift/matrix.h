#ifndef _MATRIX_H_
#define _MATRIX_H_

void identity(double **imatrix, int D);
int inverse(double **srcMatrix, int D);
int gaussjordan(double **a, double **b, int arow, int bcol);

#endif