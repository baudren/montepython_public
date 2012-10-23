cimport numpy as nm
import numpy as nm
nm.import_array()
cimport stdlib as stdlib
cimport stdio as stdio

cdef extern from "errorlist.h":
  ctypedef struct error:
    pass
    
  void stringError(char* str, error *err)
  int getErrorValue(error* err)
  int isError(error *err) 
  void purgeError(error **err) 
  void printError(void* flog,error* err)

class CError(Exception):
  def __init__(self,val,str):
    self.val=val
    self.comment=str+" "
    
  def __str__(self):

    return self.comment.strip()

cdef doError(error **err):
  cdef char estr[10000]
  if (isError(err[0])):
    stringError(estr,err[0])
    eestr = estr + " "
    er=CError(getErrorValue(err[0]),eestr)
    purgeError(err)
    
    return er
  return None

cdef extern from "wlik.h":
  ctypedef char parname[256]
  ctypedef void wlik_object
  
  void wlik_get_has_cl(wlik_object *self, int has_cl[6],error **err)
  void wlik_get_lmax(wlik_object *self, int lmax[6],error **err)
  double wlik_compute(wlik_object* self, double* cl_and_pars,error **err)
  void wlik_cleanup(wlik_object** pself)
  wlik_object* wlik_wmap_init(char *dirname, int ttmin, int ttmax, int temin, int temax, int use_gibbs, int use_lowl_pol, error **err)

cdef class wlik:
  cdef wlik_object* celf
  cdef error *_err,**err
  cdef int ndim  

  def __init__(self,char* dirname,int ttmin,int ttmax,int temin,int temax,int use_gibbs,int use_lowl_pol):
    self.celf=NULL
    self._err = NULL
    self.err = &self._err
    

    self.celf = wlik_wmap_init(dirname, ttmin, ttmax, temin, temax, use_gibbs, use_lowl_pol, self.err)
    er=doError(self.err)
    if er:
      raise er
    lmax = self.get_lmax()
    nn = nm.sum(nm.array(lmax)+1)
    self.ndim = nn 
    
        
  def __call__(self,pars):
    pars_2d = nm.atleast_2d(pars)
    if pars_2d.shape[1]!=self.ndim:
      raise Exception("Bad shape (expecting (-1,%d) got (%d,%d))"%(self.ndim,pars_2d.shape[0],pars_2d.shape[1]))
    res = nm.zeros(pars_2d.shape[0],dtype=nm.double)
    i=0
    for apars in pars_2d:
      pars_proxy=nm.PyArray_ContiguousFromAny(apars,nm.NPY_DOUBLE,1,1)
      res[i] = wlik_compute(self.celf,<double*> nm.PyArray_DATA(pars_proxy),self.err)
      er=doError(self.err)
      if er:
        raise er
      i+=1
    return res
    
  def __dealloc__(self):
    if self.celf!=NULL:
      wlik_cleanup(&(self.celf))
      
  def get_has_cl(self):
    cdef int has_cl[6]
    wlik_get_has_cl(self.celf, has_cl,self.err)
    er=doError(self.err)
    if er:
      raise er
    hcl = ""
    for i in range(6):
      hcl+=("0","1")[has_cl[i]]
    return hcl
  property has_cl:
    def __get__(self):
      return self.get_has_cl()
  
  
  def get_lmax(self):
    cdef int lmax[6]
    wlik_get_lmax(self.celf, lmax,self.err)
    er=doError(self.err)
    if er:
      raise er
    lm = ()
    for i in range(6):
      lm+=(lmax[i],)
    return lm
    
  property lmax:
    def __get__(self):
      return self.get_lmax()
  

