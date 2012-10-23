/*
 *  lklbs.c
 *  lowly_project
 *
 *  Created by Karim Benabed on 24/04/09.
 *  Copyright 2009 Institut d'Astrophysique de Paris. All rights reserved.
 *
 */
#include "lklbs.h"

int* lowly_get_ell(int *pnell, int *inell, int lmax,error **err) {
  int *outell;
  int nell,i;
  
  nell = *pnell;
  
  if (nell==0) {
    outell = malloc_err(sizeof(int)*(lmax+1),err);
    forwardError(*err,__LINE__,NULL);
    for(i=0;i<lmax+1;i++) {
      outell[i]=i;
    }
    nell = lmax+1;
  } else {
    outell = malloc_err(sizeof(int)*nell,err);
    forwardError(*err,__LINE__,NULL);
    memcpy(outell,inell,sizeof(int)*nell);
  }
  *pnell = nell;
  if (lmax!=0) {
    for(i=0;i<nell;i++) {
      testErrorRetVA(outell[i]>lmax,lowly_lbig,"ell[%d]=%d>lmax (%d)",*err,__LINE__,NULL,i,outell[i],lmax);
    }
  }
  return outell;
}



int lowly_get_offset_cl(int *has_cl, int* offset_cl,int nell) {
  int offset, i;
  offset = 0;

  for(i=0;i<6;i++) {
    if (has_cl==NULL || has_cl[i]==1) {
      offset_cl[i] = offset;
      offset += nell;
    } else {
      offset_cl[i] = -1;
    }
  }
  return offset;
}


distribution* init_fulllklbs_distribution(cmblkl** lkls,int nlkl, 
                                           bs_struct* rbs, 
                                           int *lmax, error **err) {
  distribution *dist;
  lklbs* self;

  self = init_fulllklbs(lkls,nlkl, rbs, lmax, err); 
  forwardError(*err,__LINE__,NULL);
  
  dist = init_distribution(self->ndim+self->xdim,self,&lklbs_lkl,&free_lklbs,NULL,err);
  forwardError(*err,__LINE__,NULL);
  
  return dist;
 
}

lklbs* init_fulllklbs(cmblkl** lkls,int nlkl, 
                                           bs_struct* rbs, 
                                           int *lmax, error **err) {
  lklbs* self;
  int cli,ilkl;
  int lndim;
  int ofx;
  int nrx,irx;
  
  self = malloc_err(sizeof(lklbs),err);
  forwardError(*err,__LINE__,NULL);
    
  self->rbs = rbs;
  self->ndim = rbs->ndim;

  self->lkls = malloc_err(sizeof(cmblkl*)*nlkl,err);
  forwardError(*err,__LINE__,NULL);
  
  self->ofx = malloc_err(sizeof(int)*nlkl,err);
  forwardError(*err,__LINE__,NULL);
  
  memcpy(self->lkls,lkls,sizeof(cmblkl*)*nlkl);
  self->nlkl = nlkl;
  
  lndim = 0;
  ofx = 0;
  self->xdim=0;
  nrx = 0;
  for(ilkl=0;ilkl<self->nlkl;ilkl++) {
    nrx += lkls[ilkl]->xdim;
  }
  
  if (nrx!=0) {
    self->rx = malloc_err(sizeof(int)*nrx,err);
    forwardError(*err,__LINE__,NULL);    
  } 
  self->nrx = nrx;
  
  irx = 0;
  
  for(ilkl=0;ilkl<self->nlkl;ilkl++) {
    int nd;
    int xi;
    int ipos;
      
    cmblkl_check_lmax(self->lkls[ilkl],lmax,err);
    forwardError(*err,__LINE__,NULL);
    cmblkl_check_xnames(self->lkls[ilkl],ilkl,err);
    forwardError(*err,__LINE__,NULL);
    for (xi=0;xi<lkls[ilkl]->xdim;xi++) {
      ipos = lklbs_get_par_id(self,lkls[ilkl]->xnames[xi],err);
      forwardError(*err,__LINE__,NULL);
      self->rx[irx+xi] = ipos;
    }
    self->ofx[ilkl] = irx;
    irx += lkls[ilkl]->xdim;
    nd = lkls[ilkl]->ndim + lkls[ilkl]->xdim;
    if(lndim<nd) {
      lndim=nd;
    }
    self->ofx[ilkl] = ofx;
    ofx += lkls[ilkl]->xdim;
  }    
  
  self->tot_cl = 0;
  for(cli=0;cli<6;cli++) {
    if (lmax[cli]!=-1) {
      self->offset_lmax[cli] = self->tot_cl;
      self->tot_cl += lmax[cli] + 1;
    } else {
      self->offset_lmax[cli] = -1;
    }
  }
  
  self->cl_theo=malloc_err(sizeof(double)*(lndim+self->tot_cl),err);
  forwardError(*err,__LINE__,NULL);
  self->cl_select=self->cl_theo + self->tot_cl;
  
  
  return self; 
  
}


int lklbs_get_par_id(lklbs* self,extraname name, error **err) {
  int i;
  void *old;

  if (self->rbs->xnames!=NULL) {
    // check whether the parameter is not one of the bs parameters
    for(i=0;i<self->ndim;i++) {
      if (strcmp(self->rbs->xnames[i],name)==0) {
        return i-self->ndim;
      }
    }
  }
  if (self->xdim==0) {
    self->xnames = malloc_err(sizeof(extraname),err);
    sprintf(self->xnames[0],"%s",name);
    self->xdim = 1;
    return 0;
  }
  for(i=0;i<self->xdim;i++) {
    if (strcmp(self->xnames[i],name)==0) {
      return i;
    }
  }
  old = self->xnames;
  self->xnames = malloc_err(sizeof(extraname)*(self->xdim+1),err);
  forwardError(*err,__LINE__,-1);
  memcpy(self->xnames,old,self->xdim*sizeof(extraname));
  sprintf(self->xnames[self->xdim],"%s",name);
  self->xdim++;
  free(old);
  return self->xdim-1;
}



double lklbs_lkl(void* pelf, double* pars, error **err) {
  lklbs* self;
  double res;
  double* cls;
  int cli,ilkl;
  cmblkl *llkl;
  char spars[200000];
  int i;
  
  self=pelf;
  
  /*sprintf(spars,"pars (%d,%d): ",self->ndim,self->xdim );
  for(i=0;i<self->ndim+self->xdim;i++) {
    sprintf(spars,"%s %15.10f",spars,pars[i]);
  }*/
  //_DEBUGHERE_("%s",spars);
  
  self->rbs->bs_compute(self->rbs->bs,pars,self->cl_theo,err);
  if (isError(*err)) {
    char *problem;
    int iii;
    problem = malloc(sizeof(char)*(50+10*self->ndim));
    sprintf(problem,"Problem in cl code at :\n    ");
    for (iii=0;iii<self->ndim;iii++) {
      sprintf(problem,"%s %g",problem,pars[iii]);
    }
    sprintf(problem,"%s\n",problem);
    //_DEBUGHERE_("%s",problem);
    free(problem);
  }
  forwardError(*err,__LINE__,0);
  
  res=0;
  
  for (ilkl=0;ilkl<self->nlkl;ilkl++) {
    llkl = self->lkls[ilkl];
    cls = cmblkl_select_cls(llkl,self);
    //_DEBUGHERE_("ilkl %d xdim %d ndim %d nbins %d",ilkl,llkl->xdim,llkl->ndim,llkl->nbins);
    if (llkl->xdim!=0) {
      int ii;
      for(ii=0;ii<llkl->xdim;ii++) {
        //_DEBUGHERE_("ilkl %d ixdim %d self->ndim %d ofx+i %d rx %d %g ->",ilkl,ii,self->ndim,self->ofx[ilkl]+ii,self->rx[self->ofx[ilkl]+ii],pars[self->ndim+self->rx[self->ofx[ilkl]+ii]]);
        cls[llkl->nbins+ii] = pars[self->ndim+self->rx[self->ofx[ilkl]+ii]];
      }
    }
    res += llkl->lkl_func(llkl->lkl_data,cls,err);
    forwardError(*err,__LINE__,0);
  }
  //_DEBUGHERE_("%g",res);
  return res;
}

void free_lklbs(void **pelf) {
  lklbs* self;
  int ilkl;
  cmblkl *tmp;
  
  //_DEBUGHERE_("","");
  self=*pelf;
  for (ilkl=0;ilkl<self->nlkl;ilkl++) {
    //_DEBUGHERE_("","");
    free_cmblkl(&(self->lkls[ilkl]));
  }
  //_DEBUGHERE_("","");
  free(self->lkls);
  free(self->ofx);
  //_DEBUGHERE_("","");
  free_bs_struct(&(self->rbs));
  //_DEBUGHERE_("","");
  free(self->cl_theo);
  if (self->xdim!=0) {
    free(self->xnames);
    free(self->rx);
  } 
  //_DEBUGHERE_("","");
  free(self);
  *pelf=NULL;
}

void free_cmblkl(void **pelf) {
  cmblkl* self;
  self=*pelf;
  
  //_DEBUGHERE_("","");
  if (self->lkl_free!=NULL) {
    self->lkl_free(&self->lkl_data);  
  }
  //_DEBUGHERE_("","");
  if (self->wl!=NULL) {
    free(self->wl);
  } 
  //_DEBUGHERE_("","");
  if (self->bins !=NULL) {
    free(self->pls);
  }
  //_DEBUGHERE_("","");
  if (self->xnames!=NULL) {
    free(self->xnames);
  }
  //_DEBUGHERE_("","");
  free(self);
  *pelf = NULL;
  
}

double* cmblkl_select_cls(cmblkl *llkl,lklbs* self) {
  double* cls;
  int cli,off,lff,l;
  double one;
  double *wl, *wl0;
  int inc;
  int *offset_cl,*ell;
  int ndim,xdim;
  
  
  one=1;
  if (llkl->wl==NULL) {
    wl0 = &one;
    inc = 0;
  } else {
    wl0 = llkl->wl;
    inc = 1;
  }
  offset_cl = llkl->offset_cl;
  ndim = llkl->ndim;
  xdim = llkl->xdim;
  ell = llkl->ell;
  
  if (ndim == self->tot_cl && llkl->wl==NULL && xdim==0 && llkl->unit==1) {
    cls = self->cl_theo; 
  } else {
    cls = self->cl_select; 
    for(cli=0;cli<6;cli++) {
      off = offset_cl[cli];
      lff = self->offset_lmax[cli];
      wl = wl0;
      if (off!=-1) {
        for(l=0;l<llkl->nell;l++) {
          cls[off+l] = self->cl_theo[lff+ell[l]] * (*wl) * llkl->unit;
          //_DEBUGHERE_("off %d l %d cl_t[%d]=%g cl_u[%d]=%g wl=%g",off,l,llkl->ell[l],self->cl_theo[lff+llkl->ell[l]],cls[off+l],*wl)
          wl += inc;
        }
      }
    }
  }
  if (llkl->bins!=NULL) {
    char trans;
    int npar;
    double done,dzero;
    int one;
    
    trans='T';
    npar = llkl->nbins;
    one = 1;
    done = 1.;
    dzero = 0;
    llkl->pls[0]=10;
    //_DEBUGHERE_("cls[0]=%g pls[0]=%g bins[0]=%g",cls[0],llkl->pls[0],llkl->bins[0]);
    
    dgemv(&trans, &ndim, &npar, &done, llkl->bins, &ndim, cls, &one, &dzero, llkl->pls, &one);
    //_DEBUGHERE_("cls[0]=%g pls[0]=%g ",cls[0],llkl->pls[0]);
    return llkl->pls;
  }
  return cls;
}

cmblkl *init_cmblkl_select(void* lkl_data, posterior_log_pdf_func* lkl_func, 
                    posterior_log_free *lkl_free,
                    int *lmax,
                    int xdim, error **err) {
  cmblkl_select *selec;
  cmblkl *self;
  int mlmax,i;
  int has_cl[6];
  int np;

  selec = malloc_err(sizeof(cmblkl_select),err);
  forwardError(*err,__LINE__,NULL);
  
  np = lmax[0]+1;
  mlmax = lmax[0];
  for (i=0;i<6;i++) {
    has_cl[i]=0;
    if (lmax[i]!=-1) {
      has_cl[i] = 1;  
    }
    if (lmax[i]>mlmax) {
      mlmax = lmax[i];
    }
    selec->lmax[i] = lmax[i];
    np += lmax[i] + 1;
  }
  selec->mlmax = mlmax;
  
  selec->lkl_data = lkl_data;
  selec->lkl_func = lkl_func;
  selec->lkl_free = lkl_free;
  selec->xdim = xdim;
  selec->ndim = np;
  selec->pars = malloc_err(sizeof(double)*(np+xdim),err);
  forwardError(*err,__LINE__,NULL);
  
  self = init_cmblkl(selec,select_func,select_free,0,NULL,has_cl,mlmax,1,NULL,0,NULL,0,xdim,err);
  forwardError(*err,__LINE__,NULL);
  
  return self;
}

double select_func(void* dt, double *pars,error **err) {
  cmblkl_select *selec;
  int i,cli,zero,cr;
  double res;

  selec = dt;

  zero = 0;
  cr=0;
  for(cli=0;cli<6;cli++) {
    //_DEBUGHERE_("%d :: %d %d %d,%d",cli,selec->mlmax,selec->lmax[cli],zero,cr);
    for(i=0;i<selec->lmax[cli]+1;i++) {
      selec->pars[cr] = pars[zero+i];
      cr++;
    }
    if (selec->lmax[cli]!=-1) {
      zero += selec->mlmax+1;  
    }
  }
  for(i=0;i<selec->xdim;i++) {
    //_DEBUGHERE_("%d %d %d %d %g",selec->xdim,i,cr,zero,pars[zero+i]);
    selec->pars[cr] = pars[zero+i];
    cr++;
  }

  res = selec->lkl_func(selec->lkl_data,selec->pars,err);
  forwardError(*err,__LINE__,0);

  return res;
}

void select_free(void** dt) {
  cmblkl_select *selec;

  //_DEBUGHERE_("","");
  
  selec = *dt;
  free(selec->pars);
  //_DEBUGHERE_("","");
  if(selec->lkl_free!=NULL) {
    //_DEBUGHERE_("","");
    selec->lkl_free(&(selec->lkl_data));
  }
  //_DEBUGHERE_("","");
  
  free(selec);
  *dt = NULL;
}

cmblkl *init_cmblkl(void* lkl_data, posterior_log_pdf_func* lkl_func, 
                    posterior_log_free *lkl_free,
                    int nell,int* ell,int* has_cl,int lmax,
                    double unit,
                    double *wl, int wlselect,double *bins,int nbins,int xdim, error **err) {
  cmblkl *self;
  int *_ell,_nell;
  
  self = malloc_err(sizeof(cmblkl),err);
  forwardError(*err,__LINE__,NULL);
  
  self->unit = unit;
  
  self->lkl_data = lkl_data;
  self->lkl_func = lkl_func;
  self->lkl_free = lkl_free;
  
  self->nell = nell;
  self->ell = lowly_get_ell(&self->nell, ell, lmax, err);
  forwardError(*err,__LINE__,NULL);    
  
  _ell = self->ell;
  _nell = self->nell;
  
  self->ndim = lowly_get_offset_cl(has_cl, self->offset_cl,self->nell);
  self->xdim = xdim;
  
  self->wl=NULL;
  if (wl!=NULL) {
    self->wl = malloc_err(sizeof(double)*_nell,err);
    forwardError(*err,__LINE__,NULL);    
    if (wlselect==1) {
      memcpy(self->wl,wl,sizeof(double)*_nell);
    } else {
      int il;
      for(il=0;il<_nell;il++) {
        self->wl[il] = wl[_ell[il]];
      }
    }
  }
  
  self->bins = NULL;
  self->nbins = self->ndim;  // if I have no binning matrix, nbins = ndim
  if (nbins!=0 && bins!=NULL) {
    self->nbins = nbins;
    self->pls = malloc_err(sizeof(double)*((self->ndim+1)*self->nbins+self->xdim), err);
    self->bins = self->pls + self->nbins + self->xdim;
    forwardError(*err,__LINE__,NULL);    
    memcpy(self->bins,bins,sizeof(double)*self->ndim*self->nbins);
  }
  self->xnames = NULL;
  
  return self;
}
                    
void cmblkl_set_names(cmblkl *lkl, char **names, error **err) {
  int i;
  
  if (lkl->xnames!=NULL) {
    free(lkl->xnames);
  }
  lkl->xnames = malloc_err(sizeof(extraname)*lkl->xdim,err);
  forwardError(*err,__LINE__,);
  
  for(i=0;i<lkl->xdim;i++) {
    sprintf(lkl->xnames[i],"%s",names[i]);
  }
}

void cmblkl_check_xnames(cmblkl *self,int ii,error **err) {
  int i;
  char name_tpl[50];
  
  if (self->xnames!=NULL || self->xdim==0) {
    return;
  }
  
  sprintf(name_tpl,"LKL_%d",ii);
  self->xnames = malloc_err(sizeof(extraname)*self->xdim,err);
  forwardError(*err,__LINE__,);
  
  for(i=0;i<self->xdim;i++) {
    sprintf(self->xnames[i],"%s_extra_%d",name_tpl,i);
  }
  
}

void cmblkl_check_lmax(cmblkl *lkl,int *lmax,error **err) {
  int cli;
  int _nell,*_ell,*_offset_cl;
  
  _nell = lkl->nell;
  _ell = lkl->ell;
  _offset_cl =lkl->offset_cl;
  for(cli=0;cli<6;cli++) {
    testErrorRetVA(_ell[_nell-1] > lmax[cli] && _offset_cl[cli]!=-1,lklbs_incompatible_lmax,"lmax[%d] is %d while max ell for likelihood is %d",*err,__LINE__,,cli,lmax[cli],_ell[_nell-1] > lmax[cli]);
  }
}

void cmblkl_max_lmax(cmblkl *lkl,int *lmax, error **err) {
  int cli;
  int _nell,*_ell,*_offset_cl;
  
  _nell = lkl->nell;
  _ell = lkl->ell;
  _offset_cl =lkl->offset_cl;
  
  for(cli=0;cli<6;cli++) {
    if (_offset_cl[cli]!=-1) {
      if (lmax[cli]<_ell[_nell-1]) {
        lmax[cli] = _ell[_nell-1];
      }
    }
  }
}  

zero_bs* init_zero_bs(int *lmax, error **err) {
  zero_bs *zbs;
  int i;
  
  zbs = malloc_err(sizeof(zero_bs), err);
  forwardError(*err,__LINE__,NULL);    
  zbs->ndim = 0;
  for (i=0; i<6; i++) {
    zbs->lmax[i] = lmax[i];
    //_DEBUGHERE_("zzz lmax %d %d",lmax[i],zbs->ndim);
    zbs->ndim += lmax[i] + 1;
  }
  //_DEBUGHERE_("zzz ndim %d",zbs->ndim);
  return zbs;
}

void zero_bs_compute(void* vzbs, double* prs, double* cls, error **err){
  zero_bs *zbs;
  zbs = vzbs;
  memcpy(cls,prs,sizeof(double)*zbs->ndim);
  return;
}

void free_zero_bs(void **pzbs){
  zero_bs *zbs;
  zbs = *pzbs;
  free(zbs);
  *pzbs=NULL;
  
}

bs_struct *init_bs_struct(int ndim, void* bs, compute_cl_func* bs_compute, posterior_log_free* bs_free, char **_xnames, error **err) {
  bs_struct *rbs;

  rbs = malloc_err(sizeof(bs_struct), err);
  forwardError(*err,__LINE__,NULL);

  rbs->bs = bs;
  rbs->bs_compute = bs_compute;
  rbs->bs_free = bs_free;
  rbs->ndim = ndim;
  rbs->xnames = NULL;
  if (_xnames!=NULL) {
    int i;
    rbs->xnames = malloc_err(sizeof(extraname)*ndim,err);
    forwardError(*err,__LINE__,NULL);
    for (i=0;i<ndim;i++) {
      sprintf(rbs->xnames[i],"%s",_xnames[i]);
    }
  }

  return rbs;
}

void free_bs_struct(void **prbs) {
  bs_struct *rbs;

  rbs = *prbs;
  if (rbs->bs_free!=NULL) {
    rbs->bs_free(&(rbs->bs));
  }
  if (rbs->xnames!=NULL) {
    free(rbs->xnames);
  }
  free(rbs);
  *prbs = NULL;
}

// deprecated stuff
distribution* init_lklbs_distribution(int ndim,void* lkl, posterior_log_pdf_func* lkl_func, 
                                      posterior_log_free *lkl_free, 
                                      void* bs, compute_cl_func* bs_compute, 
                                      posterior_log_free* bs_free, 
                                      int nell, int* ell, int *lmax,error **err) {
  distribution *dist;
  lklbs* self;
  
  self = init_lklbs(lkl,lkl_func,lkl_free,ndim,bs,bs_compute,bs_free,nell,ell,lmax,err);
  forwardError(*err,__LINE__,NULL);
  
  dist = init_distribution(ndim+self->xdim,self,&lklbs_lkl,&free_lklbs,NULL,err);
  forwardError(*err,__LINE__,NULL);
  
  return dist;
}

distribution* init_multilklbs_distribution(int ndim,cmblkl** lkls,int nlkl, 
                                           void* bs, compute_cl_func* bs_compute, 
                                           posterior_log_free* bs_free, 
                                           int *lmax, error **err) {
  distribution *dist;
  lklbs* self;
  
  self = init_multilklbs(lkls,nlkl,ndim,bs,bs_compute,bs_free,lmax,err);
  forwardError(*err,__LINE__,NULL);
  
  dist = init_distribution(ndim+self->xdim,self,&lklbs_lkl,&free_lklbs,NULL,err);
  forwardError(*err,__LINE__,NULL);
  
  return dist;
  
}

lklbs* init_lklbs(void* lkl, posterior_log_pdf_func* lkl_func, 
                  posterior_log_free *lkl_free, int ndim,
                  void* bs, compute_cl_func* bs_compute, 
                  posterior_log_free* bs_free, 
                  int nell, int* ell, int *lmax, error **err) {
  lklbs *self;
  int cli,clmax,mlmax,l;
  int has_cl[6];
  cmblkl *lkls;
  
  
  mlmax = -1;
  for(cli=0;cli<6;cli++) {
    if(lmax[cli]>mlmax) 
      mlmax=lmax[cli];
  }
  
  if (ell!=NULL) {
    for(l=0;l<nell;l++) {
      if (self->ell[l]-1>mlmax) {
        mlmax = self->ell[l]-1;
      }
    }
  }
  
  for(cli=0;cli<6;cli++) {
    if(lmax[cli]!=-1) {
      has_cl[cli] = 1;
      if(lmax[cli]<mlmax) 
        lmax[cli]=mlmax;
    } else{
      has_cl[cli]=-1;
    }
  }
  
  lkls = init_cmblkl(lkl, lkl_func, lkl_free, nell, ell, has_cl, mlmax,1,NULL,0,NULL,0,0, err);
  forwardError(*err,__LINE__,0);
  
  self = init_multilklbs(&lkls, 1, ndim,bs, bs_compute, bs_free, lmax, err);
  forwardError(*err,__LINE__,0);
  
  return self;
}

lklbs* init_multilklbs(cmblkl** lkls,int nlkl, 
                  int ndim,void* bs, compute_cl_func* bs_compute, 
                  posterior_log_free* bs_free, 
                  int *lmax, error **err) {
  lklbs *self;
  bs_struct *rbs;

  rbs = init_bs_struct(ndim,bs,bs_compute,bs_free,NULL,err);
  forwardError(*err,__LINE__,NULL);
  
  self =  init_fulllklbs(lkls,nlkl, rbs, lmax, err); 
  forwardError(*err,__LINE__,NULL);

  return self;
}


/*
cmblkl *beamed_distribution(cmblkl target, int xdim, void *data, beamfunc* bfunc, posterior_log_free *bfree, error **err) {
  cmblkl *dst;
  beamed *bmd;
  int has_cl[6];
  bmd = malloc_err(sizeof(beamed),err);
  forwardError(*err,__LINE__,NULL);
  
  bmd->ndim = target->ndim;
  bmd->xdim = xdim;
  bmd->target = target;
  bmd->lars = malloc_err(sizeof(double)*bmd->ndim,err);
  forwardError(*err,__LINE__,NULL);
  
  bmd->data = data;
  bmd->bfunc = bfunc;
  bmd->bfree = bfree;
  
  for(i=0;i<6;i++) {
    has_cl[i] = offset_cl[i]!=-1;
  }
  
  dst = init_cmblkl(bmd, beamed_lkl, beamed_free, bmd->target->nell, bmd->target->ell, has_cl, bmd->target->ell[bmd->target->nell-1], bmd->target->unit, bmd->target->wl, 1, bmd->target->bins, bmd->target->nbins, xdim, err);
  forwardError(*err,__LINE__,NULL);
  
  return dst;
}

double beamed_lkl(void *vbmd, double* pars, error **err) {
  beamed bmd;
  double lkl;
  
  bmd = vbmd;
  
  // apply beam
  bmd->bfunc(bmd->data, pars + bmd->ndim, pars, bmd->lars,err);
  forwardError(*err,__LINE__,0);

  // copy extra parameters
  if (bmd->target->xdim!=0) {
    memcpy(bmd->lars+bmd->ndim, pars+bmd->ndim, bmd->target->xdim*sizeof(double));
  }
  
  // compute lkl
  lkl = bmd->target->lkl_func(bmd->target->lkl_data, bmd->lars, error **err);
  forwardError(*err,__LINE__,0);

  return lkl;
}

void beamed_free(void **pbmd) {
  beamed *bmd;
  bmd = *pbmd;
  
  free_cmblkl(&(bmd->target));
  bmd->bfree(bmd->data);
  free(lars);
  free(bmd);
  *pbmd = NULL;
}
*/