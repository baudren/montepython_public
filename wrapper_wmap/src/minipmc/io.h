/* ============================================================ *
 * io.h								*
 * Martin Kilbinger, Karim Benabed 2008				*
 * ============================================================ */

#ifndef __IO_H
#define __IO_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef __PLANCK__
#include "HL2_likely/tools/errorlist.h"
#else
#include "errorlist.h"
#endif

#define io_base		-200
#define io_alloc	-1 + io_base
#define io_file		-2 + io_base
#define io_eof		-3 + io_base
#define io_null         -4 + io_base
#define io_inconsistent -5 + io_base


unsigned int numberoflines(const char *, error **err);
unsigned int numberoflines_comments(const char *name, unsigned int *ncomment, error **err);

/* The following three functions are used to read in the SNLS supernovae data */
int read_double(char** p, double *x);
int read_int(char** p, int *x);
double *readASCII(char *filename, int *fnx, int *fny, sm2_error **err);

FILE* fopen_err(const char* fname, const char *mode, error **err);
void* malloc_err(size_t sz, error **err);
void* calloc_err(size_t nl,size_t sz1, error **err);
void *realloc_err(void *ptr, size_t sz, error **err);
void* resize_err(void* orig, size_t sz_orig, size_t sz_new, int free_orig, error **err);

size_t read_line(FILE* fp, char* buff,int bmax, error **err);
void* read_any_vector(const char* fnm, size_t n, char *prs, size_t sz, error **err);
void* read_any_list(const char* fnm,size_t *n, char *prs, size_t sz, error **err);
void* read_any_list_count(const char* fnm, size_t *n, char *prs, size_t sz,
			  size_t *nlines, error **err);
double* read_double_vector(const char* fnm, size_t n, error **err);
double* read_double_list(const char* fnm, size_t *n, error **err);
int* read_int_vector(const char* fnm, size_t n, error **err);

int* read_int_list(const char* fnm, size_t *n, error **err);
float* read_float_vector(const char* fnm, size_t n, error **err);
float* read_float_list(const char* fnm, size_t *n, error **err);
long* read_long_vector(const char* fnm, size_t n, error **err);
long* read_long_list(const char* fnm,size_t *n, error **err);
void write_bin_vector(void* buf, const char* fnm, size_t n, error **err);
void* read_bin_vector(const char* fnm, size_t n, error **err);
void* read_bin_list(const char* fnm, size_t *n, error **err);

void chomp(char *x);

/* ============================================================ *
 * To calculate the program run time.				*
 * ============================================================ */
time_t start_time(FILE *FOUT);
void end_time(time_t t_start, FILE *FOUT);

#endif
