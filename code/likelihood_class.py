import os,sys

# Definition of classes for the likelihoods to inherit, for every different
# types of expected likelihoods. Included so far, models for a Clik, simple
# prior, and newdat likelihoods.

# General class that will only define the store_lkl function, common for every
# single likelihood. It will copy the content of self.path, copied from
# initialization
class likelihood():

  def __init__(self,path,command_line=False):

    self.folder = sys.path[0]+'/'
    self._read_from_file(path)
    if ((command_line is not False) and (os.path.exists(command_line.folder+'/log.dat') is False)):
      self._store_lkl_params(command_line)

  def _loglkl(self,_cosmo,data):
    raise NotImplementedError('Must implement method _loglkl() in your likelihood')

  def _store_lkl_params(self,command_line):
    log = open(command_line.folder+'log.param','a')
    tolog = open(self.path,'r')
    log.write("\n#-----Likelihood-{0}-----\n\n".format(self.__class__.__name__))
    for line in tolog:
      log.write(line)
    tolog.seek(0)
    tolog.close()
    log.close()


  def _read_from_file(self,path):
    name = self.__class__.__name__
    self.path = path
    if os.path.isfile(path):
      data_file = open(path,'r')
      for line in data_file:
	if line.find('#')==-1:
	  if line.find(name+'.')!=-1:
	    exec(line.replace(name+'.','self.'))
      data_file.seek(0)
      data_file.close()



class likelihood_prior(likelihood):

  def _loglkl(self):
    raise NotImplementedError('Must implement method __loglkl() in your likelihood')
