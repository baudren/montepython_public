import os
from likelihood_class import likelihood 

class hst(likelihood):

  def __init__(self,path='hst.data'):

    self._read_from_file(path,'hst')

    #if path is not None:
      #for line in open(path,'r'):
	#if line.find('#')==-1:
	  #if line.find('hst.')!=-1:
	    #exec(line.replace('hst.','self.'))


  def _loglkl(self,_cosmo,data):

    h   = _cosmo._h()
    loglkl = -0.5*(h-self.h)**2/(self.sigma**2)
    return loglkl
