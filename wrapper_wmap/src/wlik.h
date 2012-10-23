/*
 *  wlik.h
 *  lowly_project
 *
 *  Created by Karim Benabed on 16/03/11.
 *  Copyright 2011 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#ifndef _wlik_
#define _wlik_

#ifdef __cplusplus
extern "C" {
#endif

#include "pmc.h"
#include "lklbs.h"

#define _dealwitherr error *lerr,**err; if(_err==NULL) {lerr=NULL;err=&lerr;} else {err=_err;}

#define _forwardError(A,B,C) if(_err!=NULL) {forwardError(A,B,C);} else {quitOnError(A,B,stderr);}
#define _testErrorRetVA(A,B,C,D,E,F,...) if(_err!=NULL) {testErrorRetVA(A,B,C,D,E,F,__VA_ARGS__);} else {testErrorExitVA(A,B,C,D,E,__VA_ARGS__);}

#define _pn_size 256
typedef char parname[_pn_size];

typedef void wlik_object;
 

void wlik_get_has_cl(wlik_object *wlikid, int has_cl[6],error **err);


cmblkl* wlik_wmap_init(char *dirname, int ttmin, int ttmax, int temin, int temax, int use_gibbs, int use_lowl_pol, error **err);

void wlik_get_lmax(wlik_object *wlikid, int lmax[6],error **err);

double wlik_compute(wlik_object* wlikid, double* cl_and_pars,error **err);

// cleanup
void wlik_cleanup(wlik_object** pwlikid);

//internal
void* _wlik_dig(wlik_object* wlikid, error **err);

#ifdef __cplusplus
}
#endif

#endif