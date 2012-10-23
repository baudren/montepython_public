/*
 *  wlik.c
 *  lowly_project
 *
 *  Created by Karim Benabed on 16/03/11.
 *  Copyright 2011 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */


#include "wlik.h"

// ARE YOU STILL READING ?

// YOU HAVE BEEN WARNED !

#include <errno.h>
#include <string.h>

typedef struct {
  char tmpdir[800];
  } wmap;


void free_wmap(void **none) {
  wmap_extra_free_();
}

double wmap_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  wmap_extra_lkl_(&lkl,pars);
  return lkl;
}

cmblkl* wlik_wmap_init(char *dirname, int ttmin, int ttmax, int temin, int temax, int use_gibbs, int use_lowl_pol, error **err) {
  int bok;
  cmblkl *cing,*target;
  int mlmax;
  char pwd[5000];
  int has_cl[6],lmaxs[6];
  int *ell;
  int lmax,l,n_cl,cli;
  zero_bs* zbs;
  

  //only one wmap
  wmap_extra_only_one_(&bok);
  testErrorRet(bok!=0,-100,"wmap already initialized",*err,__LINE__,NULL);
  
  //change dir
  testErrorRetVA(getcwd(pwd,4096)==NULL,-101010,"can't get cwd name (cause = '%s')",*err,__LINE__,NULL,strerror(errno));
  testErrorRetVA(chdir(dirname)!=0,-100,"Cannot change dir to %s (cause = '%s')",*err,__LINE__,NULL,dirname,strerror(errno));

  // call wmap_init
  wmap_extra_parameter_init_(&ttmin,&ttmax,&temin,&temax,&use_gibbs,&use_lowl_pol);
  
  //change dir back
  testErrorRetVA(chdir(pwd)!=0,-100,"Cannot change dir to %s (cause = '%s')",*err,__LINE__,NULL,pwd,strerror(errno));
  
  lmax = ttmax;
  if (temax>ttmax) {
    lmax = temax;
  }
  has_cl[0]=0;
  has_cl[1]=0;
  has_cl[2]=0;
  has_cl[3]=0;
  has_cl[4]=0;
  has_cl[5]=0;
  lmaxs[0]=-1;
  lmaxs[1]=-1;
  lmaxs[2]=-1;
  lmaxs[3]=-1;
  lmaxs[4]=-1;
  lmaxs[5]=-1;

  if (ttmin<=ttmax && ttmax>0) {
    has_cl[0]=1;
  }
  if (temin<=temax && temax>0) {
    has_cl[1]=1;
    has_cl[2]=1;
    has_cl[3]=1;
  }

  ell = malloc_err(sizeof(double)*(lmax+1),err);
  forwardError(*err,__LINE__,NULL);

  for(l=0;l<=lmax;l++) {
    ell[l] = l;
  }
  cing = init_cmblkl(NULL, &wmap_lkl, 
                     &free_wmap,
                     lmax+1,ell,
                     has_cl,lmax,1,NULL,0,NULL,0,0,err);
  forwardError(*err,__LINE__,NULL);
  
  n_cl = 0;
  for(cli=0;cli<6;cli++) {
    n_cl += (lmax+1)*has_cl[cli];
  }
  
  cmblkl_max_lmax(cing,lmaxs,err);
  forwardError(*err,__LINE__,NULL);
    
  zbs = init_zero_bs(lmaxs, err);
  forwardError(*err,__LINE__,NULL);
  target = init_multilklbs_distribution(n_cl , &cing,1,
                                        zbs, &zero_bs_compute, &free_zero_bs, lmaxs, err);
  forwardError(*err,__LINE__,NULL);

  free(ell);
  
  return target;
}

void wlik_get_has_cl(wlik_object *wlikid, int has_cl[6],error **_err) {
  distribution *target;
  lklbs *lbs;
  int cli;
  _dealwitherr;

  lbs = _wlik_dig(wlikid,err);
  _forwardError(*err,__LINE__,);
  for(cli=0;cli<6;cli++) {
    //fprintf(stderr," %d %d ",cli,lbs->offset_lmax[cli]);
    if (lbs->offset_lmax[cli]!=-1) {
      has_cl[cli]=1;
    } else {
      has_cl[cli]=0;
    }
  }
}

void wlik_get_lmax(wlik_object *wlikid, int lmax[6],error **_err) {
  distribution *target;
  lklbs *lbs;
  zero_bs* zbs;
  int cli;
  _dealwitherr;
  
  lbs = _wlik_dig(wlikid,err);
  _forwardError(*err,__LINE__,);
  zbs = lbs->rbs->bs;
  
  for(cli=0;cli<6;cli++) {
    lmax[cli] = zbs->lmax[cli];
  }
}


void wlik_cleanup(wlik_object** pwlikid) {
  free_distribution(pwlikid);
}

double wlik_compute(wlik_object* wlikid, double* cl_and_pars,error **_err) {
  double res;
  _dealwitherr;
  
  res = distribution_lkl(wlikid, cl_and_pars,err);
  _forwardError(*err,__LINE__,-1);
  return res;
}

void* _wlik_dig(wlik_object* wlikid, error **err) {
  distribution *target;
  target = wlikid;
  if (target->log_pdf == &combine_lkl) { 
    // return the first wlik likelihood
    int i;
    comb_dist_data* cbd;
    cbd = target->data;
    for (i=0;i<cbd->ndist;i++) {
      if (cbd->dist[i]->log_pdf == &lklbs_lkl) {
        return cbd->dist[i]->data;
      }
    }
  }
  if (target->log_pdf==&lklbs_lkl) {
    return target->data;
  }
  testErrorRet(1==1,-111,"No wlik likelihood found",*err,__LINE__,NULL);
}

