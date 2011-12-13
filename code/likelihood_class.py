import os,sys

# Definition of classes for the likelihoods to inherit, for every different
# types of expected likelihoods. Included so far, models for a Clik, simple
# prior, and newdat likelihoods.

# General class that will only define the store_lkl function, common for every
# single likelihood. It will copy the content of self.path, copied from
# initialization
class likelihood():

  def __init__(self,path,data=None):
    raise NotImplementedError('Must implement method __init__() in your likelihood')

  def _loglkl(self,_cosmo,data):
    raise NotImplementedError('Must implement method _loglkl() in your likelihood')

  def _store_lkl_params(self,data):
    pass

  def _read_from_file(self,path,name):
    if os.path.isfile(path):
      for line in open(path,'r'):
	if line.find('#')==-1:
	  if line.find(name+'.')!=-1:
	    exec(line.replace(name+'.','self.'))



class likelihood_clik(likelihood):

  def __init__(self,path):
    raise NotImplementedError('Must implement method __init__() in your likelihood')

