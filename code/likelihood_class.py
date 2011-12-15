import os,sys
from collections import OrderedDict as od

# Definition of classes for the likelihoods to inherit, for every different
# types of expected likelihoods. Included so far, models for a Clik, simple
# prior, and newdat likelihoods.

# General class that will only define the store_lkl function, common for every
# single likelihood. It will copy the content of self.path, copied from
# initialization
class likelihood():

  def __init__(self,path,data,command_line=False):

    self.name = self.__class__.__name__
    self.folder = os.path.abspath(data.path['MontePython'])+'/../likelihoods/'+self.name+'/'
    self._read_from_file(path,data)
    if ((command_line is not False) and (os.path.exists(command_line.folder+'/log.dat') is False)):
      self._store_lkl_params(command_line)
    

  def _loglkl(self,_cosmo,data):
    raise NotImplementedError('Must implement method _loglkl() in your likelihood')

  def _store_lkl_params(self,command_line):
    log = open(command_line.folder+'log.param','a')
    tolog = open(self.path,'r')
    log.write("\n#-----Likelihood-{0}-----\n\n".format(self.name))
    for line in tolog:
      log.write(line)
    tolog.seek(0)
    tolog.close()
    log.close()

  def _read_from_file(self,path,data):
    self.path = path
    if os.path.isfile(path):
      data_file = open(path,'r')
      for line in data_file:
	if line.find('#')==-1:
	  if line.find(self.name+'.')!=-1:
	    exec(line.replace(self.name+'.','self.'))
      data_file.seek(0)
      data_file.close()

  def _need_Class_args(self,data,dictionary):
    for key,value in dictionary.iteritems():
      try :
	data.Class_args[key]
      except KeyError:
	data.Class_args[key] = ''
      if data.Class_args[key].find(value)==-1:
	data.Class_args[key] += ' '+value+' '
    pass


# Likelihood type for prior
class likelihood_prior(likelihood):

  def _loglkl(self):
    raise NotImplementedError('Must implement method __loglkl() in your likelihood')


class likelihood_newdat(likelihood):

  def _loglkl(self):
    raise NotImplementedError('Must implement method __loglkl() in your likelihood')
