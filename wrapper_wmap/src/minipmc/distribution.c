/*
 *  distribution.c
 *  ecosstat_project
 *
 *  Created by Karim Benabed on 04/11/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */

#include "pmc.h"


/******************************************************************************/
/************** Distribution interface (should be moved somewhere else) *******/
/******************************************************************************/
distribution* init_simple_distribution(int ndim,
                                void* data, 
                                posterior_log_pdf_func* log_pdf,
                                posterior_log_free* freef,
                                error **err) {
  distribution *dist;
  
  dist = init_distribution(ndim,data,log_pdf,freef,NULL,err);
  forwardError(*err,__LINE__,NULL);
  
  return dist;
}                                
                                  

distribution* init_distribution_full(int ndim,
                                     void* data, 
                                     posterior_log_pdf_func* log_pdf,
                                     posterior_log_free* freef,
                                     simulate_func *simulate,
                                     int n_ded,
                                     retrieve_ded_func* retrieve,
                                     error **err) {
  distribution *dist;
  
  testErrorRet(n_ded!=0 && retrieve==NULL,dist_undef,
               "Target invalid, expect deduced parameters, but no function provided...",
               *err,__LINE__,NULL);
  
  dist = malloc_err(sizeof(distribution),err);
  forwardError(*err,__LINE__,NULL);
  
  dist->ndim = ndim;
  dist->n_ded = n_ded;
  dist->data = data;
  dist->log_pdf = log_pdf;
  dist->free = freef;
  dist->simulate = simulate;
  dist->retrieve = retrieve;
  dist->broadcast_mpi = NULL;
  dist->ndef = 0;
  dist->def = NULL;
  dist->pars = NULL;
  dist->dlhandle = NULL;
  dist->name = NULL;
  
  dist->f_der = NULL;
  dist->d_der = NULL;
  return dist;
}

void distribution_set_broadcast(distribution* dist, mpi_exchange_func* broadcast, error **err) {
  testErrorRet(dist==NULL,-1,"Bad distribution",*err,__LINE__,);
  dist->broadcast_mpi = broadcast;
  return;
}

distribution* init_distribution(int ndim,
                                void* data, 
                                posterior_log_pdf_func* log_pdf,
                                posterior_log_free* freef,
                                simulate_func *simulate,
                                error **err) {
  distribution  *dist;
  dist = init_distribution_full(ndim, data, log_pdf, freef, simulate, 0, NULL, err);
  forwardError(*err,__LINE__,NULL);
  return dist;
}

void distribution_set_default(distribution *dist, int ndef, int* idef, double* vdef,error **err) {
  int i;
  
  testErrorRetVA(ndef >= dist->ndim + dist->ndef,dist_undef,"Too many defaults ! (expected at much %d, got %d)",*err,__LINE__,,dist->ndim+dist->ndef,ndef);

  if (dist->ndef!=0) {
    //discard old defaults
    free(dist->def);
    free(dist->pars);
    dist->ndim += dist->ndef;
    dist->ndef = 0;
  }
  
  dist->def = malloc_err(sizeof(int)*dist->ndim,err);
  forwardError(*err,__LINE__,);
  
  memset(dist->def,0,sizeof(int)*dist->ndim);
  
  dist->pars = malloc_err(sizeof(double)*dist->ndim,err);
  forwardError(*err,__LINE__,);
  
  for(i=0;i<ndef;i++) {
    testErrorRetVA(idef[i]>=dist->ndim,dist_undef,"Too many defaults ! (expected at much %d, got %d)",*err,__LINE__,,dist->ndim,idef[i]);
    testErrorRetVA(dist->def[idef[i]]==1,dist_undef,"par %d already set",*err,__LINE__,,idef[i]);
    
    dist->def[idef[i]] = 1;
    dist->pars[idef[i]] = vdef[i]; 
  }
  
  dist->ndef = ndef;
  dist->ndim -= ndef;
  
}

void distribution_set_default_name(distribution *dist, int ndef, char** namedef, double* vdef,error **err) {
  int* idef;
  
  idef = distribution_get_names(dist,ndef,namedef,0,err);
  forwardError(*err,__LINE__,);
  distribution_set_default(dist,ndef,idef,vdef,err);
  forwardError(*err,__LINE__,);
  free(idef);
}
int distribution_get_name(distribution *dist,char* name,error **err) {
  int i,j;
  
  testErrorRet(dist->name==NULL,dist_undef,"The distribution has no parameter name defined",*err,__LINE__,-1);
  j=0;
  for (i=0;i<dist->ndim + dist->ndef;i++) {
    if (dist->ndef>0 && dist->def[i]==1) {
      continue;
    } 
    if (strcmp(name,dist->name[i])==0) {
      return j;
    }
    j++;
  }

  j=-1;
  for (i=0;i<dist->n_ded;i++) {
    if (strcmp(name,dist->name[i])==0) {
      return j;
    }
    j--;
  }

  *err = addErrorVA(dist_undef, "Cannot find parameter %s", *err,__LINE__,name);
  return -1;
}

int* distribution_get_names(distribution *dist, int nnames, char** name, int includeded, error **err) {
  int* idef;
  int i;
  
  idef = malloc_err(sizeof(int)*nnames,err);
  forwardError(*err,__LINE__, NULL);
  
  for(i=0;i<nnames;i++) {
    idef[i] = distribution_get_name(dist,name[i],err);
    forwardError(*err,__LINE__, NULL);
    testErrorRetVA(idef[i]<0 && includeded==0,dist_undef,"Asking for deduced parameter (%s)",*err,__LINE__,NULL,name[i]);
  }
  
  return idef;
}


void distribution_set_names(distribution *dist,char** name, error **err) {
  int i;
  
  testErrorRet(dist->ndef!=0,dist_undef,"cannot set name after having set default parameters",*err,__LINE__,);
  if (dist->name!=NULL) {
    free(dist->name);
  }
  
  dist->name = malloc_err(sizeof(_char_name)*(dist->ndim+dist->n_ded),err);
  forwardError(*err,__LINE__,);
  
  for(i=0;i<dist->ndim+dist->n_ded;i++) {    
    sprintf(dist->name[i],"%s",name[i]);
  }
  
  return;
}

const double * distribution_fill_pars(distribution *dist, const double* pars, error **err) {
  const double *_pars;
  _pars = pars;
  if (dist->ndef!=0) {
    // i need to reset the parameters !
    int ir,ik;
    ir=0;
    for(ik=0;ik<dist->ndim+dist->ndef;ik++) {
      if (dist->def[ik]==0) {
        dist->pars[ik]=pars[ir];
        ir++;
      }
    }
    _pars = dist->pars;
  }
  return _pars;
}

double distribution_lkl(void* pdist, const double* pars, error **err) {
  distribution *dist;
  double res;
  const double *_pars;
  
  dist = pdist;
  testErrorRet(dist->log_pdf==NULL,dist_undef,"undefined log pdf function for distribution",*err,__LINE__,0);

  _pars = distribution_fill_pars(pdist, pars, err);
  forwardError(*err,__LINE__,0);
  
  res = dist->log_pdf(dist->data,_pars,err);
  forwardError(*err,__LINE__,0);
  return res;
}

void distribution_retrieve(const void* pdist, double* pars, error **err) {
  const distribution *dist;
  
  dist = pdist;
  testErrorRet(dist->retrieve==NULL && dist->n_ded!=0,dist_undef,
	       "Undefined retrieve function for distribution",*err,__LINE__,);
  if (dist->n_ded==0) {
    return;
  }
  dist->retrieve(dist->data,pars,err);
  forwardError(*err,__LINE__,);
}

void free_distribution(distribution **pdist) {
  distribution *dist;
  dist = *pdist;
  
  if (dist->name!=NULL) {
    free(dist->name);
  }
  
  if (dist->free!=NULL && dist->data!=NULL) {
    dist->free(&dist->data);
  }
  if (dist->ndef!=0) {
    free(dist->pars);
    free(dist->def);
  }
  if (dist->dlhandle!=NULL) {
    dlclose(dist->dlhandle);
  }
  free(dist);
  *pdist = NULL;
}

distribution * combine_distribution_simple_init(int ndim, error **err) {
  distribution *comb;
  
  comb = combine_distribution_init(ndim,0,err);
  forwardError(*err,__LINE__,NULL);
  
  return comb;
}

distribution * combine_distribution_init(int ndim, int nded, error **err) {
  comb_dist_data* cbd;
  distribution *self;
  int i;
  
  //_DEBUGHERE_("","");
  cbd = malloc_err(sizeof(comb_dist_data), err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");

  cbd->ndim = ndim;
  cbd->nded = nded;
  cbd->ndist = 0;
  cbd->dist = NULL;
  cbd->dim = NULL;
  cbd->ded = NULL;
  //_DEBUGHERE_("","");
  cbd->pars = malloc_err(sizeof(double)*(ndim+nded), err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  cbd->pded = malloc_err(sizeof(double)*(ndim+nded), err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  cbd->dummy = malloc_err(sizeof(double)*(ndim+nded), err);
  forwardError(*err,__LINE__,NULL);
  cbd->ndummy = ndim + nded;
  //_DEBUGHERE_("","");
  
  if (nded>0) {
    cbd->ded_from = malloc_err(sizeof(int)*(nded), err);
    forwardError(*err,__LINE__,NULL);
    //_DEBUGHERE_("","");
    for(i=0;i<nded;i++) {
      cbd->ded_from[i] = -1;
    }  
  } else {
    cbd->ded_from=NULL;
  }
  //_DEBUGHERE_("","");
  
  self = init_distribution_full(ndim, cbd, &combine_lkl, &combine_free, NULL, nded, &combine_retrieve, err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  
  return self;
}

void combine_free(void **pcbd) {
  comb_dist_data *cbd;
  int id;
  
  cbd = *pcbd;
  
  free(cbd->pars);
  free(cbd->pded);
  free(cbd->dummy);
  if (cbd->ded_from!=NULL) {
    free(cbd->ded_from);  
  }
  
  if(cbd->dist!=NULL) {
    for(id=0;id<cbd->ndist;id++) {
      free_distribution(&cbd->dist[id]);
      free(cbd->dim[id]);
      if (cbd->ded[id]!=NULL) {
        free(cbd->dim[id]);
      }
    }
    free(cbd->dist);
    free(cbd->dim);
    if(cbd->ded!=NULL) {
      free(cbd->ded);      
    }
  }
  free(cbd);
  *pcbd =NULL;
}


double combine_lkl(void *pcbd, const double* pars, error **err) {
  comb_dist_data *cbd;
  double res;
  int id;
  
  cbd = pcbd;
  res = 0;
  //_DEBUGHERE_("","");
  for(id=0;id<cbd->ndist;id++) {
    //_DEBUGHERE_("id %d",id);
    int ip;
    // select pars
    for(ip=0;ip<cbd->dist[id]->ndim;ip++) {
      int ni;
      ni = cbd->dim[id][ip];
      if (ni>=0) {
        cbd->pars[ip]=pars[ni];
      } else {
        cbd->pars[ip]=cbd->pded[-ni-1];
      }
      //_DEBUGHERE_("p %d %g",ip, cbd->pars[ip]);
    }
    //_DEBUGHERE_("r1 %g",res);
    
    // compute likelihood
    //_DEBUGHERE_("%g",cbd->pars[0]);
    res += distribution_lkl(cbd->dist[id],cbd->pars,err);
    forwardError(*err,__LINE__,0);
    //_DEBUGHERE_("r2 %g",res);
    //_DEBUGHERE_("","");
    
    // deal with deduced parameters
    if (cbd->dist[id]->n_ded!=0) {
      int ip;
      distribution_retrieve(cbd->dist[id],cbd->dummy,err);
      forwardError(*err,__LINE__,0);
      for(ip=0;ip<cbd->dist[id]->n_ded;ip++) {
        int ni;
        ni = cbd->dim[id][ip];
        if (ni>=0) {
          cbd->pded[ni] = cbd->dummy[ip];
        }
      }
    }
    //_DEBUGHERE_("","");
  }
  return res;
}
 
void combine_retrieve(const void *pcbd, double* pded, error **err) {
  const comb_dist_data *cbd;
  
  cbd = pcbd;
  memcpy(pded,cbd->pded,cbd->nded);
  return;
}

void add_to_combine_distribution_name(distribution *comb, distribution *addon, error **err) {
  int *dim_idx, *ded_idx;
  int i,j;
  
  testErrorRet(comb->name==NULL,dist_undef,"names undefined",*err,__LINE__,);
  testErrorRet(addon->name==NULL,dist_undef,"names undefined",*err,__LINE__,);
  
  dim_idx = malloc_err(sizeof(int)*addon->ndim,err);
  forwardError(*err,__LINE__,);

  j=0;
  for(i=0;i<addon->ndim+addon->ndef;i++) {
    if(addon->ndef!=0 && addon->def[i]==1) {
      continue;
    }
    dim_idx[j] = distribution_get_name(comb,addon->name[i],err);
    forwardError(*err,__LINE__,);
    j++;
  }
  
  ded_idx = malloc_err(sizeof(int)*addon->n_ded,err);
  forwardError(*err,__LINE__,);

  j=0;
  for(i=0;i<addon->n_ded;i++) {
    ded_idx[j] = distribution_get_name(comb,addon->name[i+addon->ndim+addon->ndef],err);
    forwardError(*err,__LINE__,);
    testErrorRetVA(ded_idx[j]>=0,dist_undef,"parameter %s is not deduced",*err,__LINE__,,addon->name[i+addon->ndim+addon->ndef]);
    ded_idx[j] = -ded_idx[j] -1;
    j++;
  }
  
  add_to_combine_distribution(comb,addon,dim_idx,ded_idx,err);
  forwardError(*err,__LINE__,);

  free(dim_idx);
  free(ded_idx);
}
  

void add_to_combine_distribution(distribution *comb, distribution *addon, int *dim_idx, int *ded_idx, error **err) {
  comb_dist_data *cbd;
  int i;
  
  //_DEBUGHERE_("","");
  
  testErrorRet(comb->log_pdf != &combine_lkl,dist_type,"Bad distribution", *err, __LINE__,);
  cbd = comb->data;
  //_DEBUGHERE_("","");
  
  // test whether I am compatible with the already registered distribs
  if (dim_idx!=NULL) {
    //_DEBUGHERE_("","");
    for(i=0;i<addon->ndim;i++) {
      int dix;
      dix = dim_idx[i];
      testErrorRetVA(dix>=cbd->ndim,dist_type,"addon need parameter %d, but combine_lkl only has %d dims",*err,__LINE__,,dix,cbd->ndim);
      if (dix<0) {
        testErrorRetVA(-dix-1>=cbd->nded,dist_type,"addon need parameter %d, but combine_lkl only has %d ded",*err,__LINE__,,-dix-1,cbd->nded);
        testErrorRetVA(cbd->ded_from[-dix-1]==-1,dist_type,"addon need parameter %d, but combine_lkl has no way to compute it",*err,__LINE__,,-dix-1);
      }
    }
  }
  //_DEBUGHERE_("","");
  
  if (ded_idx!=NULL) {
    //_DEBUGHERE_("","");
    for(i=0;i<addon->n_ded;i++) {
      if (ded_idx[i]>=0) {
        testErrorRetVA(cbd->ded_from[ded_idx[i]]!=-1,dist_type,"addon want to provide ded %d which is already provided by addon %d",*err,__LINE__,,ded_idx[i],cbd->ded_from[ded_idx[i]]);
      }
    }  
  }
  //_DEBUGHERE_("","");
  
  if(cbd->ndist%10 == 0) { //need resize
    //_DEBUGHERE_("","");
    cbd->dist = resize_err(cbd->dist,sizeof(distribution*)*cbd->ndist, sizeof(distribution*)*(cbd->ndist+10),1,err);
    forwardError(*err,__LINE__,);  
    cbd->dim = resize_err(cbd->dim,sizeof(int*)*cbd->ndist, sizeof(distribution*)*(cbd->ndist+10),1,err);
    forwardError(*err,__LINE__,);  
    cbd->ded = resize_err(cbd->ded,sizeof(int*)*cbd->ndist, sizeof(distribution*)*(cbd->ndist+10),1,err);
    forwardError(*err,__LINE__,);  
    //_DEBUGHERE_("","");
  }
  
  //_DEBUGHERE_("","");
  
  cbd->dist[cbd->ndist] = addon;
  //_DEBUGHERE_("","");
  
  cbd->dim[cbd->ndist] = malloc_err(sizeof(int)*addon->ndim, err);
  //_DEBUGHERE_("","");
  forwardError(*err,__LINE__,);
  if (dim_idx==NULL) {
    //_DEBUGHERE_("","");
    for(i=0;i<addon->ndim;i++) {
      cbd->dim[cbd->ndist][i]=i;
    }
    //_DEBUGHERE_("","");
  } else {
    //_DEBUGHERE_("","");
    memcpy(cbd->dim[cbd->ndist],dim_idx,sizeof(int)*addon->ndim);    
  }
  //_DEBUGHERE_("","");
  
  cbd->ded[cbd->ndist] = NULL;
  //_DEBUGHERE_("","");
  if (addon->n_ded>0) {
    //_DEBUGHERE_("","");
    cbd->ded[cbd->ndist] = malloc_err(sizeof(int)*addon->n_ded, err);
    forwardError(*err,__LINE__,);
    if (ded_idx==NULL) {
      for(i=0;i<addon->n_ded;i++) {
        cbd->ded[cbd->ndist][i]=i;
      }
    } else {
      memcpy(cbd->ded[cbd->ndist],ded_idx,sizeof(int)*addon->ndim);    
    }
  }
  //_DEBUGHERE_("","");
  
  if (addon->n_ded>cbd->ndummy) {
    free(cbd->dummy);
    cbd->ndummy = addon->n_ded;
    cbd->dummy = malloc_err(sizeof(double)*cbd->ndummy,err);
    forwardError(*err,__LINE__,);
  }
  //_DEBUGHERE_("","");
  
  for(i=0;i<addon->n_ded;i++) {
    if (ded_idx[i]>=0) {
      cbd->ded_from[ded_idx[i]] = cbd->ndist;
    }
  }
  //_DEBUGHERE_("","");
  
  cbd->ndist++;
  //_DEBUGHERE_("","");
  
}

typedef struct {
  double *mean,*tmp;
  double *std;
  int ndim;
  double logdets2;
} ezgauss;

void ezgauss_free(void **ping) {
  ezgauss * ing;
  ing = *ping;
  free(ing->std);
  free(ing->mean);
  free(ing);
  *ping = NULL;
}

double ezgauss_log_pdf(void* ping, double* pars, error **err) {
  ezgauss *ing;
  int i;
  double log_CN;
  char uplo,trans,diag;
  int xinc;
  
  ing = ping;
  log_CN=0;
  
  for(i=0;i<ing->ndim;i++) {
    ing->tmp[i] = pars[i] - ing->mean[i];
  }
  
  uplo  = 'L';
  trans = 'N';
  diag = 'N';
  xinc = 1;
  dtrsv(&uplo,&trans,&diag,&ing->ndim,ing->std,&ing->ndim,ing->tmp,&xinc);
  
  for(i=0;i<ing->ndim;i++) {
    log_CN += ing->tmp[i]*ing->tmp[i];     
  }
  return - 0.5 * (log_CN) - ing->logdets2;
}

distribution* ezgauss_init(size_t ndim, double *mean, double *Sig, error **err) {
  ezgauss *ing;
  distribution *ding;
  int plus;
  int info;
  double det;
  char uplo;
  int i;
  
  ing = malloc_err(sizeof(ezgauss), err);
  forwardError(*err,__LINE__,NULL);
  
  ing->mean = malloc_err(sizeof(double)*ndim,err);
  forwardError(*err,__LINE__,NULL);
  ing->tmp = malloc_err(sizeof(double)*ndim,err);
  forwardError(*err,__LINE__,NULL);
  ing->std = malloc_err(sizeof(double)*ndim*ndim,err);
  forwardError(*err,__LINE__,NULL);
  ing->ndim = ndim;

  memcpy(ing->mean,mean,sizeof(double)*ndim);
  memcpy(ing->std,Sig,sizeof(double)*ndim*ndim);

  uplo = 'L';
  dpotrf(&uplo,&ing->ndim,ing->std,&ing->ndim,&info);
  testErrorRetVA(info!=0,-1616165,"Could not cholesky decompose using dpotrf (%d)",*err,__LINE__,NULL,info);

  det = 1;
  for (i = 0; i < ing->ndim*ing->ndim; i+=(ing->ndim+1)) {
    det *= ing->std[i];
  }
  ing->logdets2=log(det);

  ding = init_distribution(ndim, ing, &ezgauss_log_pdf, &ezgauss_free, NULL,err);
  forwardError(*err,__LINE__,NULL);

  return ding;
}





distribution *add_gaussian_prior(distribution *orig, int ndim, int *idim, double* loc, double *var, error **err) {
  distribution *dres;
  int i;
  double *vvar;
  
  vvar = malloc_err(sizeof(double)*ndim*ndim,err);
  forwardError(*err,__LINE__,NULL);

  memset(vvar,0,sizeof(double)*ndim*ndim);
  
  for(i=0;i<ndim;i++) {
    vvar[i*ndim+i] = var[i];    
  }
  
  dres = add_gaussian_prior_2(orig,ndim,idim,loc,vvar,err);
  forwardError(*err,__LINE__,NULL);

  free(vvar);
  return dres;
}

distribution *add_gaussian_prior_2(distribution *orig, int ndim, int *idim, double* loc, double *var, error **err) {
  distribution *dres,*dmv;
  //mvdens *mv;
  
  //_DEBUGHERE_("","");
  
  //mv = mvdens_alloc(ndim,err);
  //forwardError(*err,__LINE__,NULL);
  //memcpy(mv->mean,loc,sizeof(double)*ndim);
  //memcpy(mv->std,var,sizeof(double)*ndim*ndim);
  
  //_DEBUGHERE_("","");
  
  //dmv = init_distribution(ndim, mv, &mvdens_log_pdf_void, &mvdens_free_void, NULL,err);
  //forwardError(*err,__LINE__,NULL);
  
  //_DEBUGHERE_("","");
  
  dmv = ezgauss_init(ndim,loc,var,err);
  forwardError(*err,__LINE__,NULL);
  
  if (orig->log_pdf == &combine_lkl) {
    //_DEBUGHERE_("","");
    dres = orig;
  } else {
    //_DEBUGHERE_("","");
    dres = combine_distribution_init(orig->ndim, orig->n_ded, err);
    forwardError(*err,__LINE__,NULL);
    //_DEBUGHERE_("","");
    add_to_combine_distribution(dres, orig, NULL, NULL, err);
    forwardError(*err,__LINE__,NULL);
    //_DEBUGHERE_("","");
    
  }
  //_DEBUGHERE_("","");
  
  add_to_combine_distribution(dres, dmv, idim, NULL, err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("","");
  return dres;
}

distribution *add_gaussian_prior_name(distribution *orig, int ndim, char **iname, double* loc, double *var, error **err) {
  int *idim;
  distribution *dres;
  
  idim = distribution_get_names(orig,ndim,iname,0,err);
  forwardError(*err,__LINE__,NULL);
  dres = add_gaussian_prior(orig,ndim,idim,loc,var,err);
  forwardError(*err,__LINE__,NULL);
  free(idim);
  return dres;
}

distribution *add_gaussian_prior_2_name(distribution *orig, int ndim, char **iname, double* loc, double *var, error **err) {
  int *idim;
  distribution *dres;
  
  idim = distribution_get_names(orig,ndim,iname,0,err);
  forwardError(*err,__LINE__,NULL);
  dres = add_gaussian_prior_2(orig,ndim,idim,loc,var,err);
  forwardError(*err,__LINE__,NULL);
  free(idim);
  return dres;
}

/*
void distribution_set_derivative_func(distribution *dist, first_derivative_func* f_der, second_derivative_func* d_der, error **err) {
  dist->d_der = d_der;
  dist->f_der = f_der;
}

double distribution_first_derivative(distribution *dist, int idir, double* pars, error **err) {
  double dres;
  
  
  if (dist->f_der != NULL) {
    const double *_pars;
    int _idir;
    int ik;
    _pars = distribution_fill_pars(dist, pars, err);
    forwardError(*err,__LINE__,0);

    if (dist->ndef!=0) {
       // i need to reset the parameters !
       _idir = idir;
       for(ik=0;ik<dist->ndim+dist->ndef && _idir>-1;ik++) {
         if (dist->def[ik]==0) {
           _idir--;
         }
       }
       _idir = ik-1;
    }
    dres = dist->f_der(dist->data, _idir, _pars, err);
    forwardError(*err,__LINE__,0);
  } else {
    double h;
    double errn;
    
    h = pars[idir]/1000.;
    if (h==0) {
      h = .00001; // il fadrait sans doute faire un truc plus malin ici....
    }
    dres = nd_dfridr(distribution_lkl, idir, pars, h, dist, &errn, err);
    forwardError(*err,__LINE__,0);
  }
  return dres;
}

typedef struct {
  double *pars,*pos,*along;
  int cp;
  distribution *target;
} along_struct;

double along_lkl(void* val, const double *pars, error **err) {
  along_struct *al;
  int i;
  double r;
  
  al = val;
  for(i=0;i<al->target->ndim;i++) {
    al->pars[i] = al->pos[i] + pars[0]*al->along[i];
  }
  r = distribution_lkl(al->target, al->pars, err);
  forwardError(*err,__LINE__,-1);
  return r;
}
void along_free(void** pal) {
  along_struct *al;
  
  al = *pal;
  free(al->pars);
  if(al->cp ==1) {
    free(al->pos);
  }
  free(al);
  *pal = NULL;
}

distribution *along_init(distribution *target, double *pars, double *along, int cp,error **err) {
  along_struct *al;
  distribution *dst;
  
  al = malloc_err(sizeof(along_struct),err);
  forwardError(*err,__LINE__,NULL);
  
  al->target = target;
  al->cp = cp;
  
  if (cp==1) {
    al->pos = malloc_err(sizeof(double)*target->ndim*2,err);
    forwardError(*err,__LINE__,NULL);
    al->along = al->pos + target->ndim;
    memcpy(al->pos,pars,sizeof(double)*target->ndim);
    memcpy(al->along,along,sizeof(double)*target->ndim);
  } else {
    al->pos = pars;
    al->along = along;
  }
  al->pars = malloc_err(sizeof(double)*target->ndim,err);
  forwardError(*err,__LINE__,NULL);
  dst = init_distribution_full(1,al,along_lkl, along_free,NULL,0,NULL,err);
  forwardError(*err,__LINE__,NULL);
  return dst;
}

double distribution_deriv_along(distribution *dist, double *pars, double *along, error **err) {
  distribution *alongdist;
  double r,dzero;
  
  if (dist->f_der != NULL) {
    int i;
    r = 0;
    
    for (i=0;i<dist->ndim;i++) {
      double g;
      g = distribution_first_derivative(dist,i,pars,err);
      forwardError(*err,__LINE__,-1); 
      r +=  g * along[i];
    }
    _DEBUGHERE_("deriv : %g",r);
    return r;
  }
  alongdist = along_init(dist, pars, along, 0,err);
  forwardError(*err,__LINE__,-1); 
  dzero = 0;
  r = distribution_first_derivative(alongdist,0,&dzero,err);
  forwardError(*err,__LINE__,-1); 
  free_distribution(&alongdist);
  _DEBUGHERE_("deriv : %g",r);
  return r;
} 

double distribution_second_derivative(distribution *dist, int idir,int jdir, double* pars, error **err) {
  double dres;
  
  
  if (dist->d_der != NULL) {
    const double *_pars;
    int _idir,_jdir;
    int ik;
    _pars = distribution_fill_pars(dist, pars, err);
    forwardError(*err,__LINE__,0);

    if (dist->ndef!=0) {
       // i need to reset the parameters !
       _idir = idir;
       for(ik=0;ik<dist->ndim+dist->ndef && _idir>-1;ik++) {
         if (dist->def[ik]==0) {
           _idir--;
         }
       }
       _idir = ik-1;
    }
    if (dist->ndef!=0) {
       // i need to reset the parameters !
       _jdir = jdir;
       for(ik=0;ik<dist->ndim+dist->ndef && _idir>-1;ik++) {
         if (dist->def[ik]==0) {
           _jdir--;
         }
       }
       _jdir = ik-1;
    }
    dres = dist->d_der(dist->data, _idir,_jdir, _pars, err);
    forwardError(*err,__LINE__,0);
  } else {
    double hx,hy;
    double errn;
    
    hx = pars[idir]/1000.;
    hy = pars[jdir]/1000.;
    
    if (hx==0) {
      hx = .00001; // il fadrait sans doute faire un truc plus malin ici....
    }
    if (hy==0) {
      hy = .00001; // il fadrait sans doute faire un truc plus malin ici....
    }
    dres = nd_dfridr2(distribution_lkl, idir,jdir, pars, hx,hy, dist, &errn, err);
    forwardError(*err,__LINE__,0);
  }
  return dres;
}

double* distribution_second_derivative_matrix(distribution *dist, double *pars, error **err) {
  int id, jd;
  double *hess;
  int ndim;
  double hs;

  ndim = dist->ndim;
  hess = malloc_err(sizeof(double)*dist->ndim*dist->ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  if (dist->d_der!=NULL) {
    for(id=0;id<dist->ndim;id++) {
      for(jd=id;jd<dist->ndim;jd++) {
        double hs;
        hs = distribution_second_derivative(dist,id,jd,pars,err);
        forwardError(*err,__LINE__,NULL);
        hess[id*dist->ndim+jd] = hs;
        hess[jd*dist->ndim+id] = hs;
      }
    }    
  } else {
    // first compute diagonal
    for(id=0;id<dist->ndim;id++) {
      hs = distribution_second_derivative(dist,id,id,pars,err);
      forwardError(*err,__LINE__,NULL);
      hess[id*dist->ndim+id] = hs;
    }
    for(id=0;id<dist->ndim;id++) {
      for(jd=id+1;jd<dist->ndim;jd++) {
        double hx,hy;
        double errn;

        hx = sqrt(-1/hess[id*ndim+id])/100.;
        hy = sqrt(-1/hess[jd*ndim+jd])/100.;
        
        hs = nd_dfridr2(distribution_lkl, id,jd, pars, hx,hy, dist, &errn, err);
        forwardError(*err,__LINE__,0);
        hess[id*dist->ndim+jd] = hs;
        hess[jd*dist->ndim+id] = hs;
      }
    }    
  }

  return hess;
}
*/