#include <stdio.h>
#include <stdlib.h>

#include "mem.h"

int mem4D_alloc_float(float *****arr, int octaves, int scales, int h, int w)
{
	int j;
	int mem = 0;
	(*arr) = (float ****)calloc(octaves, sizeof(float ***));
	for (j=0; j<octaves; j++)
		mem += mem3D_alloc_float((*arr)+j, scales, h, w);
	
	return mem;
}

void mem4D_free_float(float ****arr, int octaves, int scales)
{
	int i;
	if (arr)
	{
		for (i=0; i<octaves; i++)
			mem3D_free_float(arr[i], scales);
		free(arr);
		arr = NULL;
	}
}

int mem3D_alloc_float(float ****arr, int h, int w, int o)
{
	int j;
	(*arr)		= (float ***)calloc(h, sizeof(float **));
	for (j=0; j<h; j++)
		mem2D_alloc_float((*arr)+j, w, o);

	return h*w*o*sizeof(float);
}

void mem3D_free_float(float ***arr, int h)
{
	int i;
	if (arr)
	{
		for (i=0; i<h; i++)
			mem2D_free_float(arr[i]);
		free(arr);
		arr = NULL;
	}
}

int mem2D_alloc_float(float ***arr, int h, int w)
{
	int j;
	(*arr)		= (float **)calloc(h, sizeof(float *));
	(*arr)[0]	= (float *)calloc(h*w, sizeof(float));
	for (j=1; j<h; j++)
		(*arr)[j] = (*arr)[j-1] + w;

	return h*w*sizeof(float);
}

void mem2D_free_float(float **arr)
{
	if (arr)
	{
		if (arr[0])
			free(arr[0]);
		free(arr);
		arr = NULL;
	}

}

int mem2D_alloc_int(int ***arr, int h, int w)
{
	int j;
	(*arr)		= (int **)calloc(h, sizeof(int *));
	(*arr)[0]	= (int *)calloc(h*w, sizeof(int));
	for (j=1; j<h; j++)
		(*arr)[j] = (*arr)[j-1] + w;

	return h*w*sizeof(int);
}

void mem2D_free_int(int **arr)
{
	if (arr)
	{
		if (arr[0])
			free(arr[0]);
		free(arr);
		arr = NULL;
	}
}

int mem2D_alloc_double(double ***arr, int h, int w)
{
	int j;
	(*arr)		= (double **)calloc(h, sizeof(double *));
	(*arr)[0]	= (double *)calloc(h*w, sizeof(double));
	for (j=1; j<h; j++)
		(*arr)[j] = (*arr)[j-1] + w;

	return h*w*sizeof(double);
}

void mem2D_free_double(double **arr)
{
	if (arr)
	{
		if (arr[0])
			free(arr[0]);
		free(arr);
		arr= NULL;
	}
}