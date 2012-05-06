#ifndef _MEM_H_
#define _MEM_H_


int mem4D_alloc_float(float *****arr, int octaves, int scales, int h, int w);
void mem4D_free_float(float ****arr, int octaves, int scales);
int mem3D_alloc_float(float ****arr, int h, int w, int o);
void mem3D_free_float(float ***arr, int h);
int mem2D_alloc_float(float ***arr, int h, int w);
void mem2D_free_float(float **arr);
int mem2D_alloc_int(int ***arr, int h, int w);
void mem2D_free_int(int **arr);
int mem2D_alloc_double(double ***arr, int h, int w);
void mem2D_free_double(double **arr);

#endif