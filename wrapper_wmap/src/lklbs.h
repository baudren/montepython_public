/*
 *  lklbs.h
 *  lowly_project
 *
 *  Created by Karim Benabed on 24/04/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#include "pmc.h"

#ifndef __LKLBS__
#define __LKLBS__

#define TT_off 0
#define EE_off 1
#define BB_off 2
#define TE_off 3
#define TB_off 4
#define EB_off 5

#define _extra_size 256
typedef char extraname[_extra_size];

// definitions
typedef void compute_cl_func(void*, double*, double*, error **); 

typedef struct {
  void *lkl_data;
  posterior_log_pdf_func* lkl_func; 
  posterior_log_free *lkl_free;
  int lmax[6];
  int mlmax;
  int xdim,ndim;
  double *pars;
} cmblkl_select;

typedef struct {
  void *lkl_data;
  posterior_log_pdf_func* lkl_func; 
  posterior_log_free *lkl_free;
  int nell,ndim,nbins,xdim;
  int* ell;
  int offset_cl[6];
  double *wl;
  double *bins,*pls;
  double unit;
  extraname *xnames;
  } cmblkl;

typedef struct {
  int ndim;
  void* bs;
  compute_cl_func* bs_compute;
  posterior_log_free* bs_free;
  extraname *xnames;
} bs_struct;


typedef struct {
  cmblkl **lkls;
  int nlkl;
  bs_struct* rbs;
  double *cl_theo,*cl_select;
  int offset_lmax[6];
  int *ell,*ofx,*rx;
  int nell,ndim,tot_cl,nrx,xdim;
  extraname *xnames;
} lklbs;


//lklbs funcs
distribution* init_fulllklbs_distribution(cmblkl** lkls,int nlkl, 
                                           bs_struct* rbs, 
                                           int *lmax, error **err);

lklbs* init_fulllklbs(cmblkl** lkls,int nlkl, 
                                           bs_struct* rbs, 
                                           int *lmax, error **err);

void free_lklbs(void **pelf);

double lklbs_lkl(void* pelf, double* pars, error **err);


// bs support
bs_struct *init_bs_struct(int ndim, void* bs, compute_cl_func* bs_compute, posterior_log_free* bs_free, char **_xnames, error **err);
void free_bs_struct(void **prbs);

//cmblkl support
cmblkl *init_cmblkl_select(void* lkl_data, posterior_log_pdf_func* lkl_func, 
                    posterior_log_free *lkl_free,
                    int *lmax,
                    int xdim, error **err);
double select_func(void* dt, double *pars,error **err);
void select_free(void** dt);

cmblkl *init_cmblkl(void* lkl_data, posterior_log_pdf_func* lkl_func, 
                    posterior_log_free *lkl_free,
                    int nell,int* ell,int* has_cl,int lmax,double unit,double *wl,int wlselect,
                    double *bins,int nbins, int xdim,error **err);

void free_cmblkl(void **self);
void cmblkl_check_lmax(cmblkl *lkl,int *lmax,error **err);
double* cmblkl_select_cls(cmblkl *llkl,lklbs* self);
void cmblkl_max_lmax(cmblkl *lkl,int *lmax, error **err);
void cmblkl_set_names(cmblkl *lkl, char **names, error **err);
void cmblkl_check_xnames(cmblkl *self,int ii,error **err);


// support and deprecated
int lklbs_get_par_id(lklbs* self,extraname name, error **err);

lklbs* init_lklbs(void* lkl, posterior_log_pdf_func* lkl_func, 
                  posterior_log_free *lkl_free, int ndim,
                  void* bs, compute_cl_func* bs_compute, 
                  posterior_log_free* bs_free, 
                  int nell, int* ell, int *lmax,error **err);

lklbs* init_multilklbs(cmblkl** lkls,int nlkl,int ndim, 
                       void* bs, compute_cl_func* bs_compute, 
                       posterior_log_free* bs_free, 
                       int *lmax, error **err);

distribution* init_lklbs_distribution(int ndim,void* lkl, posterior_log_pdf_func* lkl_func, 
                                      posterior_log_free *lkl_free, 
                                      void* bs, compute_cl_func* bs_compute, 
                                      posterior_log_free* bs_free, 
                                      int nell, int* ell, int *lmax,error **err);
distribution* init_multilklbs_distribution(int ndim,cmblkl** lkls,int nlkl, 
                                           void* bs, compute_cl_func* bs_compute, 
                                           posterior_log_free* bs_free, 
                                           int *lmax, error **err);



// minimal bs
typedef struct {
  int ndim; 
  int lmax[6];
} zero_bs;

zero_bs* init_zero_bs(int *lmax, error **err);
void zero_bs_compute(void* zbs, double* prs, double* cls, error **err);
void free_zero_bs(void **pzbs);




/*typedef beamfunc(void* data, double* bdata, double* cls, double *bcls, error **);

typedef struct {
  int ndim;
  int xdim;
  cmblkl target;
  double *lars;
  void *data;
  beamfunc *bfunc;
  posterior_log_free *bfree;
} beamed;
*/
#define lklbs_base              -49000
#define lklbs_zero              -1 + lklbs_base
#define lklbs_incompatible_lmax -2 + lklbs_base

#define clowly_base           -9000
#define lowly_chol           -1  + clowly_base
#define lowly_unkorder       -4  + clowly_base
#define lowly_lbig           -12 + clowly_base


#endif